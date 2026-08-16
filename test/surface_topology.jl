using LinearAlgebra
using MarineHydro
using Test

"""
    structured_plate_grid(; stations, strips)

A flat structured plate, matching the fixture used by the boundary-layer tests.
Every panel carries a private copy of its vertices, exactly as the hull
importers do, so this exercises the welding path rather than bypassing it.
"""
function structured_plate_grid(; length_x = 5.0, width_y = 1.0, stations = 12,
        strips = 5)
    x = range(0.0, length_x; length = stations + 1)
    y = range(-width_y / 2, width_y / 2; length = strips + 1)
    vertex_rows = NTuple{3, Float64}[]
    face_rows = NTuple{4, Int}[]
    center_rows = NTuple{3, Float64}[]
    normal_rows = NTuple{3, Float64}[]
    areas = Float64[]
    radii = Float64[]
    strip_lists = [Int[] for _ in 1:strips]
    for j in 1:strips, i in 1:stations
        corners = ((x[i], y[j], 0.0), (x[i + 1], y[j], 0.0),
            (x[i + 1], y[j + 1], 0.0), (x[i], y[j + 1], 0.0))
        base = length(vertex_rows)
        append!(vertex_rows, corners)
        push!(face_rows, (base, base + 1, base + 2, base + 3))
        push!(center_rows, ((x[i] + x[i + 1]) / 2, (y[j] + y[j + 1]) / 2, 0.0))
        push!(normal_rows, (0.0, 0.0, 1.0))
        area = (x[i + 1] - x[i]) * (y[j + 1] - y[j])
        push!(areas, area)
        push!(radii, sqrt(area))
        push!(strip_lists[j], length(face_rows))
    end
    matrix(rows) = permutedims(reduce(hcat, [collect(row) for row in rows]))
    mesh = MarineHydro.Mesh(matrix(vertex_rows), matrix(face_rows),
        matrix(center_rows), matrix(normal_rows), areas, radii,
        length(vertex_rows), length(face_rows))
    return StructuredPanelGrid(mesh, strip_lists), stations, strips
end

kvlcc2_surface_paths() = [
    joinpath(@__DIR__, "..", "validation", "gothenburg2010", "data", "KVLCC2",
        "kvlcc_bow1.dat"),
    joinpath(@__DIR__, "..", "validation", "gothenburg2010", "data", "KVLCC2",
        "kvlcc2_stn1.dat"),
]

@testset "Surface topology" begin
    @testset "Welding recovers shared vertices" begin
        grid, stations, strips = structured_plate_grid()
        mesh = grid.mesh
        # The importer convention: no vertex is shared before welding.
        @test mesh.nvertices == 4 * mesh.nfaces
        welded, count = weld_vertices(mesh)
        # A structured plate has (stations+1) x (strips+1) distinct corners.
        @test count == (stations + 1) * (strips + 1)
        @test length(welded) == mesh.nvertices
        @test sort(unique(welded)) == 1:count
        # A positive tolerance is required.
        @test_throws ArgumentError weld_vertices(mesh; tolerance = 0.0)
    end

    @testset "Edge counts on a structured plate" begin
        grid, stations, strips = structured_plate_grid()
        topology = build_surface_topology(grid.mesh)
        interior = interior_edges(topology)
        boundary = boundary_edges(topology)
        # Interior edges: those between panels along and across the plate.
        @test length(interior) == (stations - 1) * strips + stations * (strips - 1)
        @test length(boundary) == 2 * (stations + strips)
        @test topology.nedges == length(interior) + length(boundary)
        @test all(topology.edge_cells[edge, 2] != 0 for edge in interior)
        @test all(topology.edge_cells[edge, 2] == 0 for edge in boundary)
        @test !any(==(:nonmanifold), topology.edge_kind)
        # Every panel has four real sides on a plate.
        @test all(topology.cell_edges .> 0)
        # Interior panels have four neighbours, corners have two.
        neighbour_counts = [length(topology.neighbours[panel])
                            for panel in 1:grid.mesh.nfaces]
        @test maximum(neighbour_counts) == 4
        @test minimum(neighbour_counts) == 2
    end

    @testset "Connectivity is a single component" begin
        grid, _, _ = structured_plate_grid()
        topology = build_surface_topology(grid.mesh)
        components, count = connected_components(topology)
        @test count == 1
        @test all(==(1), components)
    end

    @testset "Geometric conservation law" begin
        grid, _, _ = structured_plate_grid()
        mesh = grid.mesh
        topology = build_surface_topology(mesh)
        metrics = build_surface_metrics(mesh, topology)
        # Sum of outward edge normals weighted by edge length vanishes for a
        # closed cell. This is what catches sign and orientation errors.
        for panel in 1:mesh.nfaces
            total = zeros(2)
            for side in 1:4
                edge = topology.cell_edges[panel, side]
                edge == 0 && continue
                total .+= collect(edge_area_vector_for(topology, metrics, edge, panel))
            end
            @test norm(total) < 1e-12 * mesh.areas[panel]
        end
    end

    @testset "Local basis is orthonormal and tangent" begin
        grid, _, _ = structured_plate_grid()
        mesh = grid.mesh
        topology = build_surface_topology(mesh)
        metrics = build_surface_metrics(mesh, topology)
        for panel in 1:mesh.nfaces
            tangent = metrics.tangent[panel, :]
            binormal = metrics.binormal[panel, :]
            normal = mesh.normals[panel, :]
            @test norm(tangent)≈1 atol=1e-12
            @test norm(binormal)≈1 atol=1e-12
            @test dot(tangent, binormal)≈0 atol=1e-12
            @test dot(tangent, normal)≈0 atol=1e-12
            @test dot(binormal, normal)≈0 atol=1e-12
        end
    end

    @testset "Edge rotations are orthogonal and antisymmetric" begin
        grid, _, _ = structured_plate_grid()
        mesh = grid.mesh
        topology = build_surface_topology(mesh)
        metrics = build_surface_metrics(mesh, topology)
        for edge in interior_edges(topology)
            rotation = metrics.edge_rotation[edge, :, :]
            @test rotation' * rotation≈I atol=1e-12
            # On a planar pair the two panels see exactly opposite normals
            # once the right one is rotated into the left one's basis.
            left = topology.edge_cells[edge, 1]
            right = topology.edge_cells[edge, 2]
            left_normal = collect(edge_normal_for(topology, metrics, edge, left))
            right_normal = collect(edge_normal_for(topology, metrics, edge, right))
            @test rotation * right_normal≈-left_normal atol=1e-12
        end
        @test_throws ArgumentError edge_normal_for(topology, metrics,
            first(interior_edges(topology)), 0)
    end

    @testset "Degenerate sides are dropped" begin
        # A triangle stored as a quadrilateral with a repeated vertex must
        # contribute three edges, not four.
        vertices = [0.0 0.0 0.0; 1.0 0.0 0.0; 0.0 1.0 0.0; 0.0 1.0 0.0]
        faces = [0 1 2 3]
        centers = [1/3 1/3 0.0]
        normals = [0.0 0.0 1.0]
        mesh = MarineHydro.Mesh(vertices, faces, centers, normals, [0.5], [0.8], 4, 1)
        topology = build_surface_topology(mesh)
        @test topology.nedges == 3
        @test count(==(0), topology.cell_edges[1, :]) == 1
        @test all(!=(:interior), topology.edge_kind)
    end

    @testset "Waterline boundary classification" begin
        # A single panel with its top edge on z = 0 and the rest below.
        vertices = [0.0 0.0 0.0; 1.0 0.0 0.0; 1.0 0.0 -1.0; 0.0 0.0 -1.0]
        faces = [0 1 2 3]
        centers = [0.5 0.0 -0.5]
        normals = [0.0 1.0 0.0]
        mesh = MarineHydro.Mesh(vertices, faces, centers, normals, [1.0], [0.7], 4, 1)
        topology = build_surface_topology(mesh)
        @test count(==(:waterline), topology.edge_kind) == 1
        @test count(==(:open), topology.edge_kind) == 3
    end

    @testset "Mirrored zero is welded" begin
        # The mirror operation turns 0.0 into -0.0; those must hash together.
        grid, _, _ = structured_plate_grid(strips = 2)
        mesh = grid.mesh
        flipped = copy(mesh.vertices)
        for row in axes(flipped, 1)
            iszero(flipped[row, 2]) && (flipped[row, 2] = -0.0)
        end
        mirrored = MarineHydro.Mesh(flipped, mesh.faces, mesh.centers, mesh.normals,
            mesh.areas, mesh.radii, mesh.nvertices, mesh.nfaces)
        @test weld_vertices(mirrored)[2] == weld_vertices(mesh)[2]
    end

    @testset "Analytic hull grids weld into one surface" begin
        for grid in (surface_piercing_ellipsoid_panel_grid((3.0, 1.0, 0.8);
                longitudinal_panels = 12, girth_panels = 6),
            wigley_hull_panel_grid(2.5, 0.25, 0.156;
                longitudinal_panels = 16, vertical_panels = 4))
            topology = build_surface_topology(grid.mesh)
            _, count = connected_components(topology)
            @test count == 1
            @test !any(==(:nonmanifold), topology.edge_kind)
            # The centreplane and keel must be interior after mirroring, so the
            # only boundary left is the waterline.
            @test all(kind -> kind in (:interior, :waterline, :open),
                topology.edge_kind)
            metrics = build_surface_metrics(grid.mesh, topology)
            @test all(isfinite, metrics.edge_length)
            @test all(metrics.edge_length .> 0)
        end
    end

    @testset "KVLCC2 welds into one surface" begin
        paths = kvlcc2_surface_paths()
        if !all(isfile, paths)
            @info "skipping KVLCC2 topology test; run fetch_geometry.sh first"
        else
            for shape in ((12, 7), (16, 9))
                grid = read_gothenburg2010_panel_grid(paths; target_shape = shape)
                topology = build_surface_topology(grid.mesh)
                components, component_count = connected_components(topology)
                # THE GATE: the two point-matched blocks and the mirrored halves
                # must weld into a single connected surface, which is what makes
                # the stem ordinary interior cells instead of a marching origin.
                @test component_count == 1
                @test length(components) == grid.mesh.nfaces
                @test !any(==(:nonmanifold), topology.edge_kind)
                # The immersed double body is closed except at the free
                # surface, so the waterline must be the ONLY boundary. Any
                # `:open` edge means the bow and stern blocks failed to stitch;
                # at a tolerance of 1e-7 of the ship length a short crack
                # survives near the bilge, which is what set the default.
                @test count(==(:open), topology.edge_kind) == 0
                @test count(==(:waterline), topology.edge_kind) > 0

                metrics = build_surface_metrics(grid.mesh, topology)
                # The geometric conservation law must hold exactly, warped
                # hull panels included: this is what free-stream preservation
                # of the finite-volume scheme rests on.
                for panel in 1:grid.mesh.nfaces
                    total = zeros(2)
                    for side in 1:4
                        edge = topology.cell_edges[panel, side]
                        edge == 0 && continue
                        total .+= collect(
                            edge_area_vector_for(topology, metrics, edge, panel))
                    end
                    @test norm(total) < 1e-12 * grid.mesh.areas[panel]
                end
            end
        end
    end
end
