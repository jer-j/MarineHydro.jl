using LinearAlgebra
using MarineHydro
using Test

@testset "Cross-flow drag" begin
    grid = wigley_hull_panel_grid(4.0, 0.4, 0.25; longitudinal_panels = 16,
        vertical_panels = 4)
    mesh = grid.mesh
    all_attached = falses(mesh.nfaces)
    all_separated = trues(mesh.nfaces)

    @testset "Sectional reduction" begin
        sections = sectional_crossflow_geometry(grid, all_separated)
        @test issorted(sections.station_x)
        @test all(>(0), sections.draft)
        @test all(>(0), sections.spacing)
        @test all(==(1), sections.separated_fraction)
        # The station spacings must tile the hull exactly once.
        span = maximum(sections.station_x) - minimum(sections.station_x)
        @test sum(sections.spacing)≈span rtol=1e-12
        @test all(==(0), sectional_crossflow_geometry(grid,
            all_attached).separated_fraction)
    end

    @testset "Attached flow sheds nothing" begin
        loads = crossflow_drag_loads(grid, all_attached, 0.2, 0.0; rho = 1000.0)
        @test loads.force == 0
        @test loads.moment == 0
    end

    @testset "The loading is quadratic and opposes the drift" begin
        drag(v) = crossflow_drag_loads(grid, all_separated, v, 0.0; rho = 1000.0).force
        @test drag(0.2) < 0
        @test drag(-0.2) > 0
        @test drag(-0.2)≈-drag(0.2) rtol=1e-12
        # Quadratic: doubling the drift quadruples the force.
        @test drag(0.4)≈4 * drag(0.2) rtol=1e-12
        # Which is exactly why it cannot contribute to the linear derivative:
        # the central difference about zero is proportional to the step, so it
        # vanishes in the limit. Halving the step must halve the slope.
        slope(step) = abs((drag(step) - drag(-step)) / (2 * step))
        @test slope(1e-6)≈2 * slope(5e-7) rtol=1e-9
        @test slope(1e-7)≈slope(1e-6) / 10 rtol=1e-9
    end

    @testset "Yaw produces a moment with no net force on a symmetric hull" begin
        # A pure rotation about midship drives equal and opposite lateral flow
        # fore and aft, so the quadratic loading cancels in force and adds in
        # moment.
        centre = (maximum(mesh.centers[:, 1]) + minimum(mesh.centers[:, 1])) / 2
        loads = crossflow_drag_loads(grid, all_separated, 0.0, 0.3; rho = 1000.0,
            x_reference = centre)
        reference = crossflow_drag_loads(grid, all_separated, 0.25, 0.0; rho = 1000.0,
            x_reference = centre)
        # Not exact: the cancellation is between discrete stations, so what
        # survives is the asymmetry of the station distribution itself.
        @test abs(loads.force) < 0.01 * abs(reference.force)
        @test loads.moment < 0
    end

    @testset "Only the separated part contributes" begin
        aft = [mesh.centers[panel, 1] < 0 for panel in 1:mesh.nfaces]
        whole = crossflow_drag_loads(grid, all_separated, 0.2, 0.0; rho = 1000.0)
        partial = crossflow_drag_loads(grid, aft, 0.2, 0.0; rho = 1000.0)
        @test abs(partial.force) < abs(whole.force)
        @test abs(partial.force) > 0
        # Scaling the published coefficient scales the loading exactly.
        doubled = crossflow_drag_loads(grid, aft, 0.2, 0.0; rho = 1000.0,
            drag_coefficient = 1.6)
        @test doubled.force≈2 * partial.force rtol=1e-12
    end

    @testset "Argument validation" begin
        @test_throws ArgumentError crossflow_drag_loads(grid, all_separated, 0.2, 0.0;
            rho = 0.0)
        @test_throws ArgumentError crossflow_drag_loads(grid, all_separated, 0.2, 0.0;
            drag_coefficient = -1.0)
    end
end
