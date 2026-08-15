using ForwardDiff
using LinearAlgebra: dot
using MarineHydro
using Test

function flat_strip_grid(panel_count=8)
    panel_length = 1 / panel_count
    vertex_rows = NTuple{3,Float64}[]
    face_rows = NTuple{4,Int}[]
    centers = Matrix{Float64}(undef, panel_count, 3)
    normals = zeros(panel_count, 3)
    areas = fill(panel_length, panel_count)
    radii = fill(hypot(panel_length / 2, 0.5), panel_count)
    for panel in 1:panel_count
        x_bow = 1 - (panel - 1) * panel_length
        x_stern = x_bow - panel_length
        first_vertex = length(vertex_rows)
        append!(vertex_rows, [
            (x_bow, -0.5, 0.0),
            (x_stern, -0.5, 0.0),
            (x_stern, 0.5, 0.0),
            (x_bow, 0.5, 0.0),
        ])
        push!(face_rows, (
            first_vertex,
            first_vertex + 1,
            first_vertex + 2,
            first_vertex + 3,
        ))
        centers[panel, :] = [(x_bow + x_stern) / 2, 0.0, 0.0]
        normals[panel, 3] = 1
    end
    vertices = reduce(vcat, [collect(row)' for row in vertex_rows])
    faces = reduce(vcat, [collect(row)' for row in face_rows])
    mesh = Mesh(
        vertices,
        faces,
        centers,
        normals,
        areas,
        radii,
        size(vertices, 1),
        panel_count,
    )
    return StructuredPanelGrid(mesh, [collect(1:panel_count)])
end

@testset "Integral boundary layers" begin
    @testset "Head closure and flat plate" begin
        coordinate = collect(range(0.05, 1.0; length=20))
        edge_speed = ones(length(coordinate))
        viscosity = 1e-6
        result = @inferred solve_head_boundary_layer(
            coordinate,
            edge_speed,
            viscosity,
        )
        flat_plate_theta = 0.036 * coordinate[end] *
            (edge_speed[end] * coordinate[end] / viscosity)^(-1 / 5)
        @test all(diff(result.momentum_thickness) .> 0)
        @test result.momentum_thickness[end] ≈ flat_plate_theta rtol=0.15
        @test all(result.shape_factor .> 1.1)
        @test all(result.skin_friction_coefficient .> 0)
        @test !any(result.separated)
        @test head_kinetic_shape_factor(1.4) > 3
        @test head_entrainment_coefficient(
            head_kinetic_shape_factor(1.4),
        ) > 0

        displacement_at_speed(speed) = solve_head_boundary_layer(
            coordinate,
            fill(speed, length(coordinate)),
            viscosity,
        ).displacement_thickness[end]
        @test isfinite(ForwardDiff.derivative(displacement_at_speed, 1.0))
        @test ForwardDiff.derivative(displacement_at_speed, 1.0) < 0

        adverse_result = solve_head_boundary_layer(
            coordinate,
            collect(range(1.0, 0.4; length=length(coordinate))),
            viscosity,
        )
        @test any(adverse_result.separated)
        first_separated = something(findfirst(adverse_result.separated))
        @test all(
            iszero,
            adverse_result.skin_friction_coefficient[first_separated:end],
        )
        @test all(isfinite, adverse_result.displacement_thickness)
        @test_throws ArgumentError solve_head_boundary_layer(
            [0.0, 1.0],
            [1.0, 1.0],
            viscosity,
        )
        @test_throws ArgumentError solve_head_boundary_layer(
            [0.1, 1.0],
            [1.0, 0.0],
            viscosity,
        )
    end

    @testset "Quasi-three-dimensional surface march" begin
        grid = flat_strip_grid()
        mesh = grid.mesh
        edge_velocity = zeros(mesh.nfaces, 3)
        edge_velocity[:, 1] .= -2.0
        result = @inferred solve_quasi3d_boundary_layer(
            grid,
            edge_velocity,
            1e-6;
            rho=1000.0,
        )
        @test result.force[1] < 0
        @test abs(result.force[2]) < eps()
        @test all(result.transpiration_velocity .> 0)
        @test all(result.displacement_thickness .> 0)
        @test !any(result.separated)

        relative_velocity = @inferred body_relative_edge_velocity(
            mesh,
            2.0,
            0.1,
            0.2,
        )
        @test all(relative_velocity[:, 1] .≈ -2.0)
        @test relative_velocity[:, 2] ≈
            .-(0.1 .+ 0.2 .* mesh.centers[:, 1])
        @test all(
            abs(dot(relative_velocity[panel, :], mesh.normals[panel, :])) < 1e-14
            for panel in 1:mesh.nfaces
        )

        lateral_shear(sway_velocity) = begin
            velocity = body_relative_edge_velocity(
                mesh,
                2.0,
                sway_velocity,
                0.0,
            )
            solve_quasi3d_boundary_layer(
                grid,
                velocity,
                1e-6;
                rho=1000.0,
            ).force[2]
        end
        @test ForwardDiff.derivative(lateral_shear, 0.0) < 0
    end

    @testset "Derivative correction bookkeeping" begin
        inviscid = WangHydrodynamicDerivatives(1, 2, 3, 4, 5, 6, 7, 8)
        correction = ViscousManeuveringDerivatives(10, 20, 30, 40)
        corrected = @inferred apply_viscous_correction(inviscid, correction)
        @test corrected == WangHydrodynamicDerivatives(1, 12, 3, 24, 5, 36, 7, 48)
        prime = @inferred nondimensionalize_viscous_derivatives(
            correction,
            2.0,
            4.0;
            rho=8.0,
        )
        @test prime.Y_v == correction.Y_v / (4 * 2^2 * 4)
        @test prime.Y_r == correction.Y_r / (4 * 2^3 * 4)
        @test prime.N_v == correction.N_v / (4 * 2^3 * 4)
        @test prime.N_r == correction.N_r / (4 * 2^4 * 4)
    end

    @testset "Viscous-inviscid coupling arguments" begin
        grid = flat_strip_grid(2)
        zero_gradient = zeros(grid.mesh.nfaces, 3)
        @test_throws ArgumentError viscous_maneuvering_correction(
            grid,
            zero_gradient,
            zero_gradient,
            zero_gradient,
            1.0,
            1e-6;
            coupling_iterations=-1,
        )
        @test_throws ArgumentError viscous_maneuvering_correction(
            grid,
            zero_gradient,
            zero_gradient,
            zero_gradient,
            1.0,
            1e-6;
            coupling_relaxation=0.0,
        )
    end

    @testset "Surface streamline tracing" begin
        grid = flat_strip_grid(8)
        velocity = zeros(grid.mesh.nfaces, 3)
        velocity[:, 1] .= -1.0
        line = @inferred trace_surface_streamline(
            grid.mesh,
            velocity,
            1;
            step_size=0.04,
            max_steps=12,
            neighbor_count=2,
        )
        @test size(line.points, 2) == 3
        @test size(line.points, 1) > 2
        @test all(diff(line.points[:, 1]) .< 0)
        @test all(abs.(line.points[:, 3]) .< 1e-14)
        @test all(line.speed .≈ 1)

        velocity[4, 1] = -0.1
        stagnation = @inferred surface_stagnation_panel(grid.mesh, velocity)
        @test stagnation.panel == 4
        @test stagnation.speed ≈ 0.1
        @test_throws ArgumentError surface_stagnation_panel(
            grid.mesh,
            velocity;
            panel_mask=falses(grid.mesh.nfaces),
        )
    end
end
