using LinearAlgebra
using MarineHydro
using Test

@testset "Virtual captive tests" begin
    grid = wigley_hull_panel_grid(2.0, 0.2, 0.125; longitudinal_panels = 12,
        vertical_panels = 3)
    mesh = grid.mesh
    forward_speed = 1.0
    kinematic_viscosity = forward_speed * 2.0 / 5.0e6
    surge = solve_rigid_body_potential(mesh, :surge)

    @testset "Coefficient fitting recovers a known polynomial" begin
        # Build a matrix whose loads are an exact cubic in v' and r', then check
        # the regression returns the coefficients that generated it. This is the
        # part of the driver that has an answer independent of any hull.
        forward = 2.0
        span = 4.0
        rho = 1000.0
        truth = (Y_v = -0.021, Y_vvv = -0.35, Y_r = 0.004, Y_rrr = -0.02,
            Y_vrr = 0.11, Y_vvr = -0.07)
        force_scale = rho / 2 * span^2 * forward^2
        moment_scale = rho / 2 * span^3 * forward^2
        runs = VirtualCaptiveRun{Float64}[]
        for sway_prime in range(-0.15, 0.15; length = 7),
            yaw_prime in range(-0.4, 0.4; length = 7)

            load = truth.Y_v * sway_prime + truth.Y_vvv * sway_prime^3 +
                   truth.Y_r * yaw_prime + truth.Y_rrr * yaw_prime^3 +
                   truth.Y_vrr * sway_prime * yaw_prime^2 +
                   truth.Y_vvr * sway_prime^2 * yaw_prime
            push!(runs,
                VirtualCaptiveRun(sway_prime * forward, yaw_prime * forward / span,
                    load * force_scale, 2 * load * moment_scale, 0.0, 0.0, 0, true))
        end
        straight = VirtualCaptiveRun(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0, true)
        fitted = fit_captive_coefficients(runs, straight, forward, span; rho)
        for name in keys(truth)
            @test getproperty(fitted.sway, name)≈getproperty(truth, name) rtol=1e-10
        end
        # The yaw column was built as twice the sway one.
        @test fitted.yaw.N_v≈2 * truth.Y_v rtol=1e-10
        @test fitted.yaw.N_vvv≈2 * truth.Y_vvv rtol=1e-10
        # An exactly representable matrix must leave no residual.
        @test fitted.residual.sway < 1e-12
        @test fitted.residual.yaw < 1e-12
    end

    @testset "A pure drift sweep fits the sway family alone" begin
        runs = [VirtualCaptiveRun(0.1 * step, 0.0, 3.0 * step, -1.0 * step, 0.0, 0.0,
                    0, true) for step in -3:3]
        straight = VirtualCaptiveRun(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0, true)
        fitted = fit_captive_coefficients(runs, straight, 1.0, 2.0)
        @test keys(fitted.sway) == (:Y_v, :Y_vvv)
        @test keys(fitted.yaw) == (:N_v, :N_vvv)
        @test_throws ArgumentError fit_captive_coefficients(
            [VirtualCaptiveRun(0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 0, true)], straight,
            1.0, 2.0)
    end

    @testset "Sweeping a hull is antisymmetric in drift" begin
        # A hull symmetric about its centreplane must give a sway force odd in
        # the drift angle, whatever the closure does in between. This catches a
        # sign or reference-point error in the assembly of the loads.
        angles = [-0.08, -0.04, 0.0, 0.04, 0.08]
        test = virtual_captive_test(grid, surge.potential_gradient, nothing, nothing,
            forward_speed, kinematic_viscosity; drift_angles = angles,
            yaw_rates = [0.0], rho = 1000.0)
        @test length(test.runs) == length(angles)
        forces = [run.sway_force - test.straight.sway_force for run in test.runs]
        moments = [run.yaw_moment - test.straight.yaw_moment for run in test.runs]
        scale = maximum(abs, forces)
        @test scale > 0
        @test forces[1]≈-forces[end] rtol=1e-6
        @test forces[2]≈-forces[end - 1] rtol=1e-6
        @test abs(forces[3]) < 1e-8 * scale
        @test moments[1]≈-moments[end] rtol=1e-6
        # Drift must produce a side force opposing the drift.
        @test test.coefficients.sway.Y_v < 0
        @test all(run -> run.converged, test.runs)
    end

    @testset "Argument validation" begin
        @test_throws ArgumentError virtual_captive_test(grid,
            surge.potential_gradient, nothing, nothing, -1.0, kinematic_viscosity)
        @test_throws ArgumentError virtual_captive_test(grid,
            surge.potential_gradient, nothing, nothing, forward_speed, -1.0)
        @test_throws ArgumentError virtual_captive_test(grid,
            surge.potential_gradient, nothing, nothing, forward_speed,
            kinematic_viscosity; drift_angles = Float64[])
    end
end
