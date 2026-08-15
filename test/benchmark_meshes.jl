using MarineHydro
using ForwardDiff
using Test

@testset "Gothenburg structured hull surfaces" begin
    mktemp() do path, io
        write(
            io,
            """TITLE = \"synthetic half hull\"
VARIABLES = \"X\", \"Y\", \"Z\"
ZONE I=2 J=2 F=POINT
-0.5 0.25 -0.5
 0.5 0.25 -0.5
-0.5 0.25  0.0
 0.5 0.25  0.0
""",
        )
        flush(io)

        patches = @inferred read_tecplot_structured_surfaces(path)
        @test length(patches) == 1
        @test size(only(patches)) == (2, 2, 3)
        @test only(patches)[1, 1, :] == [-0.5, 0.25, -0.5]
        @test only(patches)[2, 2, :] == [0.5, 0.25, 0.0]

        mesh = @inferred read_gothenburg2010_mesh(path; length_scale=2.0)
        @test mesh.nfaces == 2
        @test mesh.nvertices == 8
        @test mesh_surface_area(mesh) ≈ 4.0
        @test sort(mesh.centers[:, 2]) ≈ [-0.5, 0.5]
        @test sort(mesh.normals[:, 2]) ≈ [-1.0, 1.0]
        @test all(mesh.centers[:, 1] .≈ 0.0)
        @test all(mesh.centers[:, 3] .≈ -0.5)
        @test all(mesh.areas .≈ 2.0)

        grid = @inferred read_gothenburg2010_panel_grid(
            path;
            length_scale=2.0,
        )
        @test grid.mesh.nfaces == mesh.nfaces
        @test length(grid.strips) == 2
        @test sort(reduce(vcat, grid.strips)) == collect(1:mesh.nfaces)
        @test all(
            issorted(grid.mesh.centers[strip, 1]; rev=true)
            for strip in grid.strips
        )

        half_mesh = read_gothenburg2010_mesh(
            path;
            length_scale=2.0,
            mirror=false,
        )
        @test half_mesh.nfaces == 1
        @test only(half_mesh.normals[:, 2]) ≈ -1.0
    end

    mktemp() do path, io
        write(io, "ZONE I=2 J=2 F=POINT\n0 0 0\n")
        flush(io)
        @test_throws ArgumentError read_tecplot_structured_surfaces(path)
    end

    @testset "Maximum-section cutoff" begin
        mktemp() do path, io
            write(
                io,
                """ZONE I=3 J=3 F=POINT
-0.2 0.0 -1.0
 0.0 0.0 -1.0
 0.2 0.0 -1.0
-0.2 1.0 -1.0
 0.0 2.0 -1.0
 0.2 1.9 -1.0
-0.2 1.0  0.0
 0.0 2.0  0.0
 0.2 1.9  0.0
""",
            )
            flush(io)

            patch = only(read_tecplot_structured_surfaces(path))
            curve = @inferred structured_sectional_area_curve(patch)
            half_curve = @inferred structured_sectional_area_curve(
                patch;
                mirror=false,
            )
            @test curve.x ≈ [-0.2, 0.0, 0.2]
            @test curve.area ≈ [2.0, 4.0, 3.8]
            @test half_curve.area ≈ curve.area ./ 2

            exact = @inferred gothenburg_maximum_section(
                path;
                relative_area_tolerance=0.0,
            )
            @test exact.x ≈ 0.0
            @test exact.maximum_x ≈ 0.0
            @test exact.maximum_area ≈ 4.0
            @test exact.candidate_count == 1

            aftmost = @inferred gothenburg_maximum_section(
                path;
                length_scale=2.0,
                relative_area_tolerance=0.1,
            )
            @test aftmost.x ≈ -0.4
            @test aftmost.area ≈ 15.2
            @test aftmost.maximum_area ≈ 16.0
            @test aftmost.candidate_count == 2

            @test_throws ArgumentError gothenburg_maximum_section(
                path;
                relative_area_tolerance=1.0,
            )
        end
    end
end

@testset "Public benchmark surface formats" begin
    mktemp() do path, io
        write(
            io,
            """1
2 2 1
0 1 0 1
0.25 0.25 0.25 0.25
-1 -1 1 1
""",
        )
        flush(io)
        patches = @inferred read_plot3d_structured_surfaces(path)
        @test length(patches) == 1
        @test size(only(patches)) == (2, 2, 3)
        @test only(patches)[:, 1, 3] == [-1.0, -1.0]

        mesh = structured_surface_waterline_mesh(
            patches;
            mirror=true,
            interior_point=[0.5, 0.0, -0.5],
        )
        @test mesh.nfaces == 2
        @test maximum(mesh.vertices[:, 3]) == 0.0
        @test mesh_surface_area(mesh) ≈ 2.0
        @test sort(mesh.normals[:, 2]) ≈ [-1.0, 1.0]
    end

    mktemp() do path, io
        write(
            io,
            """solid clipped
facet normal 0 -1 0
outer loop
vertex 0 0 -1
vertex 1 0 -1
vertex 1 0 1
endloop
endfacet
facet normal 0 -1 0
outer loop
vertex 0 0 -1
vertex 1 0 1
vertex 0 0 1
endloop
endfacet
endsolid clipped
""",
        )
        flush(io)
        mesh = read_ascii_stl_waterline_mesh(
            path;
            clustering=(4, 1, 4),
            interior_point=[0.5, 0.1, -0.5],
        )
        @test mesh.nfaces >= 1
        @test maximum(mesh.vertices[:, 3]) == 0.0
        @test mesh_surface_area(mesh) ≈ 1.0
    end
end

@testset "Wigley hull benchmark" begin
    hull_length = 2.5
    breadth = 0.25
    draft = 0.156
    grid = @inferred wigley_hull_panel_grid(
        hull_length,
        breadth,
        draft;
        longitudinal_panels=40,
        vertical_panels=10,
    )
    @test grid.mesh.nfaces == 800
    @test length(grid.strips) == 20
    @test all(length(strip) == 40 for strip in grid.strips)
    @test all(
        issorted(grid.mesh.centers[strip, 1]; rev=true)
        for strip in grid.strips
    )
    exact_volume = 4hull_length * breadth * draft / 9
    @test mesh_signed_volume(grid.mesh) ≈ exact_volume rtol=4e-3
    @test maximum(grid.mesh.vertices[:, 3]) == 0.0

    volume_sensitivity = ForwardDiff.derivative(length_value -> begin
        mesh = wigley_hull_mesh(
            length_value,
            breadth,
            draft;
            longitudinal_panels=12,
            vertical_panels=4,
        )
        mesh_signed_volume(mesh)
    end, hull_length)
    @test volume_sensitivity ≈ 4breadth * draft / 9 rtol=0.01
end

@testset "Analytical ellipsoid benchmark" begin
    sphere_coefficients = @inferred ellipsoid_potential_coefficients(
        (1.0, 1.0, 1.0),
    )
    @test collect(sphere_coefficients) ≈ fill(2 / 3, 3) rtol=2e-13
    @test ellipsoid_potential_coefficients((1, 1, 1)) == sphere_coefficients

    semi_axes = (3.0, 1.0, 0.8)
    coefficients = @inferred ellipsoid_potential_coefficients(semi_axes)
    @test sum(coefficients) ≈ 2.0 rtol=2e-13
    @test coefficients.x < coefficients.y < coefficients.z

    density = 1000.0
    sphere_added_mass = @inferred ellipsoid_added_mass(
        (1.0, 1.0, 1.0);
        density,
    )
    expected_sphere_mass = density * 2pi / 3
    @test collect(sphere_added_mass) ≈ fill(expected_sphere_mass, 3) rtol=2e-13

    added_mass_sensitivity = ForwardDiff.derivative(
        axis -> ellipsoid_added_mass(
            (axis, semi_axes[2], semi_axes[3]);
            density,
        ).x,
        semi_axes[1],
    )
    @test isfinite(added_mass_sensitivity)

    volume_sensitivity = ForwardDiff.derivative(semi_major_axis -> begin
        differentiated_mesh = surface_piercing_ellipsoid_mesh(
            (semi_major_axis, semi_axes[2], semi_axes[3]);
            longitudinal_panels=8,
            girth_panels=4,
        )
        mesh_signed_volume(differentiated_mesh)
    end, semi_axes[1])
    @test isfinite(volume_sensitivity)
    @test volume_sensitivity > 0

    coarse = @inferred surface_piercing_ellipsoid_panel_grid(
        semi_axes;
        longitudinal_panels=12,
        girth_panels=6,
    )
    fine = surface_piercing_ellipsoid_panel_grid(
        semi_axes;
        longitudinal_panels=32,
        girth_panels=16,
    )
    @test coarse.mesh.nfaces == 72
    @test fine.mesh.nfaces == 512
    @test length(coarse.strips) == 6
    @test all(length(strip) == 12 for strip in coarse.strips)
    @test all(
        issorted(coarse.mesh.centers[strip, 1]; rev=true)
        for strip in coarse.strips
    )
    @test maximum(coarse.mesh.vertices[:, 3]) <= 10eps(Float64)
    @test minimum(vec(sum(
        coarse.mesh.centers .* coarse.mesh.normals;
        dims=2,
    ))) > 0

    exact_half_volume = 2pi * prod(semi_axes) / 3
    coarse_error = abs(mesh_signed_volume(coarse.mesh) / exact_half_volume - 1)
    fine_error = abs(mesh_signed_volume(fine.mesh) / exact_half_volume - 1)
    @test fine_error < coarse_error
    @test fine_error < 0.01

    @test_throws DimensionMismatch ellipsoid_potential_coefficients(
        (1.0, 1.0),
    )
    @test_throws ArgumentError ellipsoid_potential_coefficients(
        (1.0, -1.0, 1.0),
    )
    @test_throws ArgumentError ellipsoid_potential_coefficients(
        (1.0, 1.0, 1.0);
        quadrature_order=4,
    )
    @test_throws ArgumentError surface_piercing_ellipsoid_mesh(
        semi_axes;
        longitudinal_panels=3,
    )
end
