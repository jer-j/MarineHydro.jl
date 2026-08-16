using LinearAlgebra
using MarineHydro
using Test

"""
    square_plate_mesh(; side, cells)

A flat square plate as a panel mesh, with private vertices per panel exactly as
the hull importers produce. Square rather than long so that a flow crossing it
at a large angle still has a useful run length.
"""
function square_plate_mesh(; side = 5.0, cells = 40)
    nodes = range(0.0, side; length = cells + 1)
    vertex_rows = NTuple{3, Float64}[]
    face_rows = NTuple{4, Int}[]
    center_rows = NTuple{3, Float64}[]
    normal_rows = NTuple{3, Float64}[]
    areas = Float64[]
    radii = Float64[]
    for j in 1:cells, i in 1:cells
        corners = ((nodes[i], nodes[j], 0.0), (nodes[i + 1], nodes[j], 0.0),
            (nodes[i + 1], nodes[j + 1], 0.0), (nodes[i], nodes[j + 1], 0.0))
        base = length(vertex_rows)
        append!(vertex_rows, corners)
        push!(face_rows, (base, base + 1, base + 2, base + 3))
        push!(center_rows, ((nodes[i] + nodes[i + 1]) / 2,
            (nodes[j] + nodes[j + 1]) / 2, 0.0))
        push!(normal_rows, (0.0, 0.0, 1.0))
        area = (nodes[i + 1] - nodes[i]) * (nodes[j + 1] - nodes[j])
        push!(areas, area)
        push!(radii, sqrt(area))
    end
    matrix(rows) = permutedims(reduce(hcat, [collect(row) for row in rows]))
    return MarineHydro.Mesh(matrix(vertex_rows), matrix(face_rows),
        matrix(center_rows), matrix(normal_rows), areas, radii,
        length(vertex_rows), length(face_rows))
end

"""
    plate_run_length(x, y, angle, side)

Distance travelled along the flow from whichever inflow edge of the square feeds
the point, for a uniform stream at `angle` to the x axis.
"""
function plate_run_length(x, y, angle, side)
    along_x = cos(angle) > 1e-12 ? x / cos(angle) : Inf
    along_y = sin(angle) > 1e-12 ? y / sin(angle) : Inf
    return min(along_x, along_y)
end

"""
    plate_seam_offset(x, y, angle, side)

Signed distance, as a fraction of `side`, from the line where the two inflow
edges feed equal run lengths.

`plate_run_length` takes a minimum of two competing origins, so on that line it
has a corner and the single-origin flat-plate law stops describing the flow:
two layers with different histories merge there. The band around it is excluded
from the comparison because the oracle is undefined in it, not because the
solver misbehaves. At zero incidence there is only one origin and the offset is
infinite everywhere, so nothing is excluded.
"""
function plate_seam_offset(x, y, angle, side)
    (cos(angle) <= 1e-12 || sin(angle) <= 1e-12) && return Inf
    return (x / cos(angle) - y / sin(angle)) / side
end

@testset "Surface boundary layer" begin
    viscosity = 1e-6
    speed = 1.0
    side = 5.0

    @testset "Oblique flat plate at $(round(Int, rad2deg(angle))) degrees" for angle in
                                                                              (0.0,
        deg2rad(30), deg2rad(60))
        mesh = square_plate_mesh(; side, cells = 40)
        edge_velocity = zeros(mesh.nfaces, 3)
        edge_velocity[:, 1] .= speed * cos(angle)
        edge_velocity[:, 2] .= speed * sin(angle)

        layer = solve_surface_boundary_layer(mesh, edge_velocity, viscosity;
            rho = 1000.0)

        # A uniform stream over a flat plate cannot generate crossflow at any
        # angle to the mesh lines. This is the invariance the strip march did
        # not have.
        @test maximum(abs, layer.crossflow_angle) < 1e-8
        @test !any(layer.separated)

        # Compare against the turbulent flat-plate laws, measuring run length
        # along the flow. Sample well away from both inflow edges so the
        # starting transient and the corner are excluded.
        run_lengths = [plate_run_length(mesh.centers[panel, 1],
                           mesh.centers[panel, 2], angle, side)
                       for panel in 1:mesh.nfaces]
        seam_offsets = [plate_seam_offset(mesh.centers[panel, 1],
                            mesh.centers[panel, 2], angle, side)
                        for panel in 1:mesh.nfaces]
        sampled = findall((run_lengths .> 0.45 * side) .&
                          (abs.(seam_offsets) .> 0.1))
        @test length(sampled) > 50
        for panel in sampled
            reynolds = speed * run_lengths[panel] / viscosity
            expected_thickness = 0.036 * run_lengths[panel] * reynolds^(-0.2)
            expected_friction = 0.0592 * reynolds^(-0.2)
            @test layer.momentum_thickness[panel]≈expected_thickness rtol=0.10
            @test layer.skin_friction_coefficient[panel]≈expected_friction rtol=0.10
        end

        # Wall traction opposes nothing here — it acts along the imposed stream,
        # so it must be parallel to it.
        traction_angle = atan(layer.force[2], layer.force[1])
        @test traction_angle≈angle atol=1e-6
        @test layer.force[3]≈0 atol=1e-10
    end

    @testset "Rigid-motion invariance" begin
        # Rotating the whole problem in space must not change the answer. This
        # is what catches a momentum flux rotated as R*Theta instead of
        # R*Theta*R'.
        mesh = square_plate_mesh(; side, cells = 24)
        angle = deg2rad(40)
        edge_velocity = zeros(mesh.nfaces, 3)
        edge_velocity[:, 1] .= speed * cos(angle)
        edge_velocity[:, 2] .= speed * sin(angle)
        reference = solve_surface_boundary_layer(mesh, edge_velocity, viscosity;
            rho = 1000.0)

        rotation = let a = 0.7, b = -0.4
            rx = [1 0 0; 0 cos(a) -sin(a); 0 sin(a) cos(a)]
            rz = [cos(b) -sin(b) 0; sin(b) cos(b) 0; 0 0 1]
            rz * rx
        end
        turned = MarineHydro.Mesh(mesh.vertices * rotation', mesh.faces,
            mesh.centers * rotation', mesh.normals * rotation', mesh.areas,
            mesh.radii, mesh.nvertices, mesh.nfaces)
        turned_layer = solve_surface_boundary_layer(turned,
            edge_velocity * rotation', viscosity; rho = 1000.0)

        @test turned_layer.momentum_thickness≈reference.momentum_thickness rtol=1e-9
        @test turned_layer.shape_factor≈reference.shape_factor rtol=1e-9
        @test turned_layer.skin_friction_coefficient≈
              reference.skin_friction_coefficient rtol=1e-9
        @test norm(rotation * reference.force - turned_layer.force) <
              1e-8 * norm(reference.force)
    end

    @testset "Mesh-line invariance" begin
        # The same physical problem on a coarser and a finer plate must agree
        # where they overlap, and neither may depend on the mesh orientation.
        angle = deg2rad(55)
        results = map((30, 45)) do cells
            mesh = square_plate_mesh(; side, cells)
            edge_velocity = zeros(mesh.nfaces, 3)
            edge_velocity[:, 1] .= speed * cos(angle)
            edge_velocity[:, 2] .= speed * sin(angle)
            layer = solve_surface_boundary_layer(mesh, edge_velocity, viscosity;
                rho = 1000.0)
            runs = [plate_run_length(mesh.centers[panel, 1], mesh.centers[panel, 2],
                        angle, side) for panel in 1:mesh.nfaces]
            seams = [plate_seam_offset(mesh.centers[panel, 1], mesh.centers[panel, 2],
                         angle, side) for panel in 1:mesh.nfaces]
            interior = findall((runs .> 0.5 * side) .& (abs.(seams) .> 0.1))
            ratios = [layer.momentum_thickness[panel] /
                      (0.036 * runs[panel] * (speed * runs[panel] / viscosity)^(-0.2))
                      for panel in interior]
            sum(ratios) / length(ratios)
        end
        @test all(ratio -> isapprox(ratio, 1.0; rtol = 0.10), results)
        @test results[1]≈results[2] rtol=0.05
    end

    @testset "Argument validation" begin
        mesh = square_plate_mesh(; cells = 6)
        edge_velocity = zeros(mesh.nfaces, 3)
        edge_velocity[:, 1] .= 1.0
        @test_throws DimensionMismatch solve_surface_boundary_layer(mesh,
            zeros(mesh.nfaces, 2), 1e-6)
        @test_throws ArgumentError solve_surface_boundary_layer(mesh, edge_velocity,
            -1.0)
        @test_throws ArgumentError solve_surface_boundary_layer(mesh, edge_velocity,
            1e-6; rho = 0.0)
    end
end
