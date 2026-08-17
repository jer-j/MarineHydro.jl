using LinearAlgebra
using MarineHydro
using Test

@testset "Shed vorticity" begin
    grid = wigley_hull_panel_grid(4.0, 0.4, 0.25; longitudinal_panels = 20,
        vertical_panels = 5)
    mesh = grid.mesh
    sections = sectional_crossflow_geometry(grid, falses(mesh.nfaces))
    stern = minimum(sections.station_x)
    bow = maximum(sections.station_x)
    forward_speed = 2.0
    rho = 1000.0

    @testset "A body whose draft tapers recovers d'Alembert" begin
        # The momentum balance telescopes to zero only when the sectional added
        # mass vanishes at the tail, and m22 goes as the square of the DRAFT.
        # A Wigley hull tapers in beam but not in draft, so it does not qualify;
        # a surface-piercing ellipsoid does. Getting this backwards is easy and
        # would make the model look like it satisfies d'Alembert when it is
        # really just reporting a small tail section.
        tapered = surface_piercing_ellipsoid_panel_grid((3.0, 1.0, 0.8);
            longitudinal_panels = 24, girth_panels = 8)
        taper_sections = sectional_crossflow_geometry(tapered,
            falses(tapered.mesh.nfaces))
        taper_stern = minimum(taper_sections.station_x)
        closed = shed_vorticity_derivatives(taper_sections, taper_stern,
            forward_speed; rho)
        scale = rho * forward_speed * maximum(taper_sections.draft)^2
        @test abs(closed.Y_v) < 0.15 * scale
        @test closed.N_v < 0
    end

    @testset "Separation opposes the drift, at the depth it separates" begin
        for fraction in (0.05, 0.1, 0.2, 0.4)
            station = stern + fraction * (bow - stern)
            result = shed_vorticity_derivatives(sections, station, forward_speed; rho)
            # A sway velocity to starboard must be resisted.
            @test result.Y_v < 0
            # The magnitude is set by the section depth where the flow leaves,
            # so on a hull of constant draft it barely depends on where that is.
            # That insensitivity is a feature: the linear derivative does not
            # hinge on pinning the separation station precisely.
            @test result.Y_v≈-forward_speed * rho * π * maximum(sections.draft)^2 rtol=0.1
        end
    end

    @testset "Magnitude matches the low-aspect-ratio estimate" begin
        # Once separation reaches full draft the classical slender-body result
        # is Y_v' = -2 pi (T/L)^2. Check the model reproduces it.
        draft = maximum(sections.draft)
        length_scale = bow - stern
        station = stern + 0.5 * length_scale
        result = shed_vorticity_derivatives(sections, station, forward_speed; rho)
        nondimensional = result.Y_v / (rho / 2 * length_scale^2 * forward_speed)
        expected = -2π * (draft / length_scale)^2
        @test nondimensional≈expected rtol=0.25
        @test nondimensional < 0
    end

    @testset "Moving the separation station is smooth" begin
        # The mesh derivative passes through this lookup, so a small move of the
        # station must give a proportionally small change rather than a jump
        # between stations.
        base = stern + 0.3 * (bow - stern)
        step = 1e-4 * (bow - stern)
        values = [shed_vorticity_derivatives(sections, base + offset * step,
                      forward_speed; rho).Y_v for offset in -2:2]
        differences = diff(values)
        @test all(!iszero, differences)
        @test maximum(abs, differences)≈minimum(abs, differences) rtol=0.2
    end

    @testset "Argument validation" begin
        @test_throws ArgumentError shed_vorticity_derivatives(sections, stern, -1.0)
        @test_throws ArgumentError shed_vorticity_derivatives(sections,
            bow + 1.0, forward_speed)
        @test_throws ArgumentError sectional_added_mass(sections; rho = 0.0)
    end
end
