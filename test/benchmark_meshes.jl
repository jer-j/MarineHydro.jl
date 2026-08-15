using MarineHydro
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
