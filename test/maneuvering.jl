using LinearAlgebra
using MarineHydro
using Test
using ForwardDiff

"""
    wetted_box_mesh(; half_length=1.0, half_beam=0.5, draft=1.0)

Construct the five-panel wetted surface of a surface-piercing rectangular box.
The top of the box lies on ``z = 0`` and is omitted because the reflected
Rankine kernel enforces the rigid-free-surface condition.
"""
function wetted_box_mesh(; half_length=1.0, half_beam=0.5, draft=1.0)
    a = half_length
    b = half_beam
    d = draft
    vertices = [
        -a -b -d
        -a  b -d
         a  b -d
         a -b -d
        -a -b  0
        -a  b  0
         a  b  0
         a -b  0
    ]
    faces = [
        0 1 2 3
        3 2 6 7
        1 0 4 5
        1 5 6 2
        3 7 4 0
    ]
    centers = [
         0  0 -d
         a  0 -d / 2
        -a  0 -d / 2
         0  b -d / 2
         0 -b -d / 2
    ]
    normals = [
         0  0 -1
         1  0  0
        -1  0  0
         0  1  0
         0 -1  0
    ]
    areas = [4a * b, 2b * d, 2b * d, 2a * d, 2a * d]
    radii = [
        hypot(a, b),
        hypot(b, d / 2),
        hypot(b, d / 2),
        hypot(a, d / 2),
        hypot(a, d / 2),
    ]
    return Mesh(
        Float64.(vertices),
        Int.(faces),
        Float64.(centers),
        Float64.(normals),
        Float64.(areas),
        Float64.(radii),
        size(vertices, 1),
        size(faces, 1),
    )
end

@testset "Wang maneuvering derivatives" begin
    mesh = wetted_box_mesh()

    @testset "Boundary condition and stern mask" begin
        @test @inferred(wang_sway_boundary_condition(mesh)) == mesh.normals[:, 2]
        @test @inferred(maneuvering_boundary_condition(mesh, :surge)) ==
            mesh.normals[:, 1]
        @test @inferred(maneuvering_boundary_condition(mesh, :sway)) == mesh.normals[:, 2]
        x_reference = 0.25
        y_reference = -0.1
        expected_yaw =
            -(mesh.centers[:, 2] .- y_reference) .* mesh.normals[:, 1] .+
            (mesh.centers[:, 1] .- x_reference) .* mesh.normals[:, 2]
        @test maneuvering_boundary_condition(
            mesh,
            :yaw;
            x_reference,
            y_reference,
        ) == expected_yaw
        @test_throws ArgumentError maneuvering_boundary_condition(mesh, :roll)
        @test @inferred(wang_stern_mask(mesh, 0.0)) ==
            BitVector(mesh.centers[:, 1] .>= 0)
        @test wang_stern_mask(
            mesh,
            0.0;
            bow_direction=:negative_x,
        ) == BitVector(mesh.centers[:, 1] .<= 0)
        @test_throws ArgumentError wang_stern_mask(mesh, 0.0; bow_direction=:invalid)
    end

    @testset "Wang restricted-water elevation" begin
        coordinate = [-0.5, 0.0, 0.5]
        section_area = [0.0, 0.02, 0.0]
        waterline_beam = [0.0, 0.2, 0.0]
        unrestricted = @inferred wang_restricted_water_elevation(
            coordinate,
            section_area,
            waterline_beam,
            1.0,
        )
        @test unrestricted.local_speed == ones(3)
        @test unrestricted.elevation == zeros(3)
        @test iszero(unrestricted.mean_elevation)

        restricted = @inferred wang_restricted_water_elevation(
            coordinate,
            section_area,
            waterline_beam,
            1.0;
            channel_area=1.0,
            channel_surface_width=1.0,
        )
        @test restricted.local_speed[2] > 1
        @test restricted.elevation[2] < 0
        continuity_residual = restricted.local_speed[2] * (
            1.0 - section_area[2] +
            restricted.elevation[2] * (1.0 - waterline_beam[2])
        ) - 1.0
        @test abs(continuity_residual) < 1e-12
        @test restricted.mean_elevation < 0
        @test_throws ArgumentError wang_restricted_water_elevation(
            coordinate,
            section_area,
            waterline_beam,
            1.0;
            channel_area=1.0,
        )
    end

    @testset "Dual-valued mesh compatibility" begin
        function yaw_mode_from_center(center_x)
            scale = one(center_x)
            dual_mesh = Mesh(
                mesh.vertices .* scale,
                mesh.faces,
                mesh.centers .* scale,
                mesh.normals .* scale,
                mesh.areas .* scale,
                mesh.radii .* scale,
                mesh.nvertices,
                mesh.nfaces,
            )
            dual_mesh.centers[4, 1] = center_x
            return maneuvering_boundary_condition(dual_mesh, :yaw)[4]
        end

        @test ForwardDiff.derivative(yaw_mode_from_center, 0.3) ≈ 1.0
    end

    @testset "Independent sway and yaw surface integration" begin
        sway_potential = [-2.0, -1.0, 3.0, 4.0, -5.0]
        yaw_potential = [1.5, -2.5, 0.75, -1.25, 3.0]
        sway_phi_x = [1.0, -2.0, 0.5, -1.0, 2.0]
        yaw_phi_x = [-0.5, 1.25, -2.0, 0.25, 1.5]
        speed = 2.5
        rho = 1000.0
        x_reference = 0.2
        y_reference = -0.1
        mask = BitVector([true, true, false, true, false])
        result = potential_flow_hydrodynamic_derivatives(
            mesh,
            sway_potential,
            yaw_potential,
            sway_phi_x,
            yaw_phi_x,
            speed;
            rho,
            x_reference,
            y_reference,
            velocity_mask=mask,
        )

        g_v = maneuvering_boundary_condition(mesh, :sway)
        g_r = maneuvering_boundary_condition(
            mesh,
            :yaw;
            x_reference,
            y_reference,
        )
        sway_weight = g_v .* mesh.areas
        yaw_weight = g_r .* mesh.areas
        @test result.Y_vdot == rho * sum(sway_potential .* sway_weight)
        @test result.Y_rdot == rho * sum(yaw_potential .* sway_weight)
        @test result.N_vdot == rho * sum(sway_potential .* yaw_weight)
        @test result.N_rdot == rho * sum(yaw_potential .* yaw_weight)
        @test result.Y_v == -rho * speed * sum(sway_phi_x .* sway_weight .* mask)
        @test result.Y_r == -rho * speed * sum(yaw_phi_x .* sway_weight .* mask)
        @test result.N_v == -rho * speed * sum(sway_phi_x .* yaw_weight .* mask)
        @test result.N_r == -rho * speed * sum(yaw_phi_x .* yaw_weight .* mask)
    end

    @testset "Surface integration" begin
        potential = [-2.0, -1.0, 3.0, 4.0, -5.0]
        phi_x = [1.0, -2.0, 0.5, -1.0, 2.0]
        speed = 2.5
        rho = 1000.0
        mask = BitVector([true, true, false, true, false])
        result = wang_hydrodynamic_derivatives(
            mesh,
            potential,
            phi_x,
            speed;
            rho,
            velocity_mask=mask,
        )

        x = mesh.centers[:, 1]
        weight = mesh.normals[:, 2] .* mesh.areas
        velocity_weight = weight .* mask
        @test result.Y_vdot == rho * sum(potential .* weight)
        @test result.Y_rdot == rho * sum(x .* potential .* weight)
        @test result.N_vdot == result.Y_rdot
        @test result.N_rdot == rho * sum(x .^ 2 .* potential .* weight)
        @test result.Y_v == -rho * speed * sum(phi_x .* velocity_weight)
        @test result.Y_r == -rho * speed * sum(
            (potential .+ x .* phi_x) .* velocity_weight,
        )
        @test result.N_v == -rho * speed * sum(x .* phi_x .* velocity_weight)
        @test result.N_r == -rho * speed * sum(
            x .* (potential .+ x .* phi_x) .* velocity_weight,
        )

        result_twice_speed = wang_hydrodynamic_derivatives(
            mesh,
            potential,
            phi_x,
            2 * speed;
            rho,
            velocity_mask=mask,
        )
        @test result_twice_speed.Y_vdot == result.Y_vdot
        @test result_twice_speed.Y_v == 2 * result.Y_v
        @test result_twice_speed.Y_r == 2 * result.Y_r
        @test result_twice_speed.N_v == 2 * result.N_v
        @test result_twice_speed.N_r == 2 * result.N_r
    end

    @testset "Wang nondimensionalization" begin
        dimensional = WangHydrodynamicDerivatives(
            4.0,
            8.0,
            16.0,
            32.0,
            64.0,
            128.0,
            256.0,
            512.0,
        )
        L = 2.0
        U = 4.0
        rho = 8.0
        nondimensional = nondimensionalize_wang_derivatives(
            dimensional,
            L,
            U;
            rho,
        )
        @test nondimensionalize_maneuvering_derivatives(
            dimensional,
            L,
            U;
            rho,
        ) == nondimensional
        half_rho = rho / 2
        @test nondimensional.Y_vdot == dimensional.Y_vdot / (half_rho * L^3)
        @test nondimensional.Y_v == dimensional.Y_v / (half_rho * L^2 * U)
        @test nondimensional.Y_rdot == dimensional.Y_rdot / (half_rho * L^4)
        @test nondimensional.Y_r == dimensional.Y_r / (half_rho * L^3 * U)
        @test nondimensional.N_vdot == dimensional.N_vdot / (half_rho * L^4)
        @test nondimensional.N_v == dimensional.N_v / (half_rho * L^3 * U)
        @test nondimensional.N_rdot == dimensional.N_rdot / (half_rho * L^5)
        @test nondimensional.N_r == dimensional.N_r / (half_rho * L^4 * U)
        @test_throws ArgumentError nondimensionalize_wang_derivatives(
            dimensional,
            L,
            0.0;
            rho,
        )
    end

    @testset "MMG-to-Wang velocity normalization" begin
        mmg = @inferred MMGLinearHullDerivatives(-0.315, 0.083, -0.137, -0.049)
        converted = @inferred mmg_to_wang_velocity_derivatives(mmg, 20.8, 320.0)
        draft_ratio = 20.8 / 320.0
        @test converted.Y_v == draft_ratio * mmg.Y_v
        @test converted.Y_r == draft_ratio * mmg.Y_r
        @test converted.N_v == draft_ratio * mmg.N_v
        @test converted.N_r == draft_ratio * mmg.N_r
        @test ForwardDiff.derivative(
            draft -> mmg_to_wang_velocity_derivatives(mmg, draft, 320.0).Y_v,
            20.8,
        ) ≈ mmg.Y_v / 320.0
        @test_throws ArgumentError mmg_to_wang_velocity_derivatives(mmg, 0.0, 320.0)
        @test_throws ArgumentError mmg_to_wang_velocity_derivatives(mmg, 20.8, 0.0)
    end

    @testset "Indirect potential gradient" begin
        green_functions = (Rankine(), RankineReflected())
        sources = ComplexF64[0.2, -0.3, 0.4, -0.5, 0.6]
        gradient = evaluate_indirect_potential_gradient(
            green_functions,
            mesh,
            sources,
            0.0,
        )

        for axis in 1:3
            direction = zeros(3)
            direction[axis] = 1
            _, K_direction = assemble_matrices(
                green_functions,
                mesh,
                0.0;
                direct=false,
                all_normals=direction,
            )
            @test gradient[:, axis] ≈ K_direction * sources atol=1e-12 rtol=1e-12
        end
    end

    @testset "End-to-end zero-frequency solve" begin
        result = solve_wang_maneuvering(mesh, 2.0; rho=1000.0)
        normal_gradient = vec(sum(result.potential_gradient .* mesh.normals; dims=2))

        @test result.boundary_residual < 1e-10
        @test normal_gradient ≈ wang_sway_boundary_condition(mesh) atol=1e-10
        @test result.derivatives.Y_rdot == result.derivatives.N_vdot
        @test all(result.velocity_mask)
        @test all(isfinite, result.potential)
        @test all(isfinite, result.potential_gradient)

        dY_r_dU = ForwardDiff.derivative(
            speed -> solve_wang_maneuvering(mesh, speed; rho=1000.0).derivatives.Y_r,
            2.0,
        )
        @test dY_r_dU ≈ result.derivatives.Y_r / 2.0
    end

    @testset "Unit-surge potential" begin
        result = solve_rigid_body_potential(mesh, :surge)
        normal_gradient = vec(sum(
            result.potential_gradient .* mesh.normals;
            dims=2,
        ))
        @test result.boundary_residual < 1e-10
        @test normal_gradient ≈ maneuvering_boundary_condition(
            mesh,
            :surge,
        ) atol=1e-10
        @test all(isfinite, result.potential)
        @test all(isfinite, result.potential_gradient)
    end


    @testset "Independent unit-yaw zero-frequency solve" begin
        x_reference = 0.25
        result = solve_potential_flow_maneuvering(
            mesh,
            2.0;
            rho=1000.0,
            x_reference,
        )
        sway_normal_gradient = vec(sum(
            result.potential_gradients.sway .* mesh.normals;
            dims=2,
        ))
        yaw_normal_gradient = vec(sum(
            result.potential_gradients.yaw .* mesh.normals;
            dims=2,
        ))

        @test result.boundary_residuals.sway < 1e-10
        @test result.boundary_residuals.yaw < 1e-10
        @test sway_normal_gradient ≈ maneuvering_boundary_condition(
            mesh,
            :sway,
        ) atol=1e-10
        @test yaw_normal_gradient ≈ maneuvering_boundary_condition(
            mesh,
            :yaw;
            x_reference,
        ) atol=1e-10
        @test result.derivatives.Y_rdot ≈ result.derivatives.N_vdot atol=1e-10
        @test result.potentials.yaw ≈ -x_reference .* result.potentials.sway
        @test all(isfinite, result.potential_gradients.yaw)

        dN_r_dU = ForwardDiff.derivative(
            speed -> solve_potential_flow_maneuvering(
                mesh,
                speed;
                rho=1000.0,
                x_reference,
            ).derivatives.N_r,
            2.0,
        )
        @test dN_r_dU ≈ result.derivatives.N_r / 2.0
    end
end
