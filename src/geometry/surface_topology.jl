using LinearAlgebra: cross, dot, norm

"""
    row3(matrix, index)

Row `index` of an `n × 3` array as a stack-allocated `SVector`.

Subtracting two `@view`s of a matrix row builds a heap-allocated vector every
time, and these differences sit in the innermost loops of the surface solver —
gradient stencils, edge upwinding, curvature. Reading the row into an `SVector`
instead keeps the arithmetic on the stack. On KVLCC2 at 480 panels this removed
about a hundred thousand allocations per cache build.
"""
@inline row3(matrix, index::Integer) = SVector{3}(matrix[index, 1], matrix[index, 2],
    matrix[index, 3])

raw"""
    SurfaceTopology

Edge connectivity of a panel mesh, for finite-volume discretisation over the
hull surface.

Every mesh produced by the structured importers gives each panel a *private*
copy of its four vertices, so `nvertices == 4 * nfaces` and adjacency cannot be
read from `Mesh.faces`. It is instead recovered by welding coincident vertices
within a tolerance; see [`weld_vertices`](@ref).

The topology is integer-valued and is built from the primal (derivative-free)
coordinates. That is deliberate rather than a compromise: connectivity is
piecewise constant in the geometry, so its derivative vanishes almost
everywhere. Everything that carries a derivative lives in
[`SurfaceMetrics`](@ref).

# Fields

- `welded`: mesh vertex index (1-based) to welded vertex identifier.
- `edge_vertices`: `(nedges, 2)` welded identifiers, ascending.
- `edge_cells`: `(nedges, 2)` panel indices. The second entry is `0` on a
  boundary edge.
- `edge_local_side`: `(nedges, 2)` which of the owning panel's four sides this
  edge is, `0` where there is no second panel.
- `cell_edges`: `(nfaces, 4)` edge index per panel side, `0` for a collapsed
  side. Triangles are stored as quadrilaterals with a repeated vertex, so one
  side of such a panel has zero length and no edge.
- `edge_kind`: `:interior`, `:waterline`, `:open` or `:nonmanifold`.
- `neighbours` / `neighbour_cells`: compressed adjacency used by the
  least-squares gradient stencil.
"""
struct SurfaceTopology
    welded::Vector{Int}
    edge_vertices::Matrix{Int}
    edge_cells::Matrix{Int}
    edge_local_side::Matrix{Int}
    cell_edges::Matrix{Int}
    edge_kind::Vector{Symbol}
    neighbours::Vector{UnitRange{Int}}
    neighbour_cells::Vector{Int}
    nedges::Int
end

# Strip any dual/tracked number down to a plain value so that the welding
# lattice is computed on coordinates alone.
_primal(x::Real) = x
_primal(x::ForwardDiff.Dual) = _primal(ForwardDiff.value(x))

# The mirror operation negates the transverse coordinate, which turns 0.0 into
# -0.0. Those hash to different lattice cells, so normalise the sign of zero.
_unsigned_zero(x) = iszero(x) ? zero(x) : x

struct _UnionFind
    parent::Vector{Int}
end

_UnionFind(n::Integer) = _UnionFind(collect(1:n))

function _find(sets::_UnionFind, index::Integer)
    root = index
    while sets.parent[root] != root
        root = sets.parent[root]
    end
    while sets.parent[index] != root
        sets.parent[index], index = root, sets.parent[index]
    end
    return root
end

function _union!(sets::_UnionFind, left::Integer, right::Integer)
    left_root, right_root = _find(sets, left), _find(sets, right)
    left_root == right_root && return left_root
    sets.parent[right_root] = left_root
    return left_root
end

raw"""
    weld_vertices(mesh; tolerance = nothing)

Merge coincident mesh vertices and return `(welded, count)`, where `welded[i]`
is the identifier of mesh vertex `i` and `count` is the number of distinct
welded vertices.

Vertices are bucketed on a lattice of pitch `tolerance` and merged against the
twenty-seven neighbouring buckets, so points separated by less than the
tolerance are welded regardless of which bucket they land in.

The default tolerance is `1e-6` of the largest model dimension. That is set by
the multi-block hull importers rather than by round-off: on KVLCC2 the bow and
stern blocks are point-matched in the source data, but resampling leaves seam
vertices disagreeing at a few times `1e-7` of the ship length, which at `1e-7`
leaves a short crack near the bilge. The welded edge count is flat from `1e-6`
through `1e-3` of the model span, so the default sits in a wide basin — well
above the seam mismatch and three orders below the panel size.
"""
function weld_vertices(mesh::Mesh; tolerance = nothing)
    coordinates = [_unsigned_zero(_primal(mesh.vertices[vertex, axis]))
                   for vertex in 1:mesh.nvertices, axis in 1:3]
    span = maximum(
        maximum(@view(coordinates[:, axis])) - minimum(@view(coordinates[:, axis]))
        for axis in 1:3)
    pitch = isnothing(tolerance) ? 1e-6 * max(span, one(span)) : tolerance
    pitch > 0 || throw(ArgumentError("weld tolerance must be positive"))

    buckets = Dict{NTuple{3, Int}, Vector{Int}}()
    cell_of(vertex) = ntuple(
        axis -> floor(Int, coordinates[vertex, axis] / pitch), 3)
    for vertex in 1:mesh.nvertices
        push!(get!(buckets, cell_of(vertex), Int[]), vertex)
    end

    sets = _UnionFind(mesh.nvertices)
    squared_tolerance = pitch^2
    for vertex in 1:mesh.nvertices
        base = cell_of(vertex)
        for di in -1:1, dj in -1:1, dk in -1:1
            neighbours = get(buckets, (base[1] + di, base[2] + dj, base[3] + dk),
                nothing)
            isnothing(neighbours) && continue
            for other in neighbours
                other >= vertex && continue
                distance = sum(
                    (coordinates[vertex, axis] - coordinates[other, axis])^2
                    for axis in 1:3)
                distance <= squared_tolerance && _union!(sets, vertex, other)
            end
        end
    end

    identifiers = Dict{Int, Int}()
    welded = Vector{Int}(undef, mesh.nvertices)
    for vertex in 1:mesh.nvertices
        root = _find(sets, vertex)
        welded[vertex] = get!(identifiers, root, length(identifiers) + 1)
    end
    return welded, length(identifiers)
end

raw"""
    build_surface_topology(mesh; kwargs...)

Build the [`SurfaceTopology`](@ref) of a panel mesh.

# Keywords

- `weld_tolerance`: passed to [`weld_vertices`](@ref).
- `waterline`: elevation of the free surface, default `0`.
- `waterline_tolerance`: a boundary edge whose two endpoints both lie within
  this distance of `waterline` is classified `:waterline`. Default is `1e-6` of
  the model's vertical extent.
- `warn_nonmanifold`: emit a warning when an edge is shared by more than two
  panels, which the clustered STL importer can produce.

Edges shared by exactly two panels are `:interior`. Everything else is a
boundary: `:waterline` where the free surface cuts the hull, `:open` elsewhere
(a transom, or an unstitched multi-block seam), and `:nonmanifold` where more
than two panels meet. Inflow and outflow are *not* classified here — that
depends on the flow direction and is decided per solve from the sign of the
edge-normal velocity.
"""
function build_surface_topology(mesh::Mesh; weld_tolerance = nothing,
        waterline = 0.0, waterline_tolerance = nothing, warn_nonmanifold = true)
    welded, welded_count = weld_vertices(mesh; tolerance = weld_tolerance)
    heights = [_primal(mesh.vertices[vertex, 3]) for vertex in 1:mesh.nvertices]
    vertical_span = maximum(heights) - minimum(heights)
    height_tolerance = isnothing(waterline_tolerance) ?
                       1e-6 * max(vertical_span, one(vertical_span)) :
                       waterline_tolerance
    # One representative height per welded vertex; they are coincident to within
    # the weld tolerance so any contributing mesh vertex will do.
    welded_height = zeros(eltype(heights), welded_count)
    for vertex in 1:mesh.nvertices
        welded_height[welded[vertex]] = heights[vertex]
    end

    # Collect the four sides of every panel, dropping any side that collapses
    # to a point. A triangle is stored as a quadrilateral with a repeated
    # vertex, so exactly one of its sides is degenerate.
    sides = Dict{Tuple{Int, Int}, Vector{Tuple{Int, Int}}}()
    cell_edges = zeros(Int, mesh.nfaces, 4)
    for panel in 1:mesh.nfaces
        corners = ntuple(k -> welded[mesh.faces[panel, k] + 1], 4)
        for side in 1:4
            first_corner = corners[side]
            second_corner = corners[mod1(side + 1, 4)]
            first_corner == second_corner && continue
            key = minmax(first_corner, second_corner)
            push!(get!(sides, key, Tuple{Int, Int}[]), (panel, side))
        end
    end

    nedges = length(sides)
    edge_vertices = zeros(Int, nedges, 2)
    edge_cells = zeros(Int, nedges, 2)
    edge_local_side = zeros(Int, nedges, 2)
    edge_kind = Vector{Symbol}(undef, nedges)
    nonmanifold_count = 0

    for (edge, (key, owners)) in enumerate(sides)
        edge_vertices[edge, 1], edge_vertices[edge, 2] = key
        for (slot, (panel, side)) in enumerate(owners)
            slot > 2 && break
            edge_cells[edge, slot] = panel
            edge_local_side[edge, slot] = side
            cell_edges[panel, side] = edge
        end
        if length(owners) == 2
            edge_kind[edge] = :interior
        elseif length(owners) > 2
            edge_kind[edge] = :nonmanifold
            nonmanifold_count += 1
            # Sides beyond the first two are recorded on the panel but the edge
            # is treated as a boundary, so no flux crosses it.
            for (panel, side) in owners[3:end]
                cell_edges[panel, side] = edge
            end
            edge_cells[edge, 2] = 0
            edge_local_side[edge, 2] = 0
        else
            on_waterline = all(
                abs(welded_height[identifier] - waterline) <= height_tolerance
                for identifier in key)
            edge_kind[edge] = on_waterline ? :waterline : :open
        end
    end

    if warn_nonmanifold && nonmanifold_count > 0
        @warn "surface topology: $(nonmanifold_count) edges are shared by more " *
              "than two panels and are treated as boundaries" nonmanifold_count
    end

    neighbours, neighbour_cells = _cell_adjacency(mesh.nfaces, edge_cells, edge_kind)
    return SurfaceTopology(welded, edge_vertices, edge_cells, edge_local_side,
        cell_edges, edge_kind, neighbours, neighbour_cells, nedges)
end

function _cell_adjacency(nfaces, edge_cells, edge_kind)
    lists = [Int[] for _ in 1:nfaces]
    for edge in axes(edge_cells, 1)
        edge_kind[edge] === :interior || continue
        left, right = edge_cells[edge, 1], edge_cells[edge, 2]
        push!(lists[left], right)
        push!(lists[right], left)
    end
    neighbour_cells = Int[]
    neighbours = Vector{UnitRange{Int}}(undef, nfaces)
    for panel in 1:nfaces
        start = length(neighbour_cells) + 1
        append!(neighbour_cells, lists[panel])
        neighbours[panel] = start:length(neighbour_cells)
    end
    return neighbours, neighbour_cells
end

"""
    boundary_edges(topology)
    interior_edges(topology)

Indices of the boundary and interior edges of a [`SurfaceTopology`](@ref).
"""
boundary_edges(topology::SurfaceTopology) =
    findall(!=(:interior), topology.edge_kind)
interior_edges(topology::SurfaceTopology) =
    findall(==(:interior), topology.edge_kind)

raw"""
    connected_components(topology)

Partition panels into groups connected through interior edges. A hull that has
welded correctly returns exactly one component; more than one means a
multi-block seam failed to stitch.
"""
function connected_components(topology::SurfaceTopology)
    nfaces = size(topology.cell_edges, 1)
    sets = _UnionFind(nfaces)
    for edge in 1:topology.nedges
        topology.edge_kind[edge] === :interior || continue
        _union!(sets, topology.edge_cells[edge, 1], topology.edge_cells[edge, 2])
    end
    labels = Dict{Int, Int}()
    components = Vector{Int}(undef, nfaces)
    for panel in 1:nfaces
        root = _find(sets, panel)
        components[panel] = get!(labels, root, length(labels) + 1)
    end
    return components, length(labels)
end

raw"""
    SurfaceMetrics

Geometric quantities of a panel mesh needed by a finite-volume surface solver,
in the element type of the mesh so that they carry derivatives.

`tangent` and `binormal` are the local surface basis ``\hat e_1``, ``\hat e_2``
of each panel.

`edge_area_vector[e, slot, :]` is the outward edge normal *times* the edge
length, for the panel in `edge_cells[e, slot]`, expressed in that panel's own
basis. Storing the product rather than a unit normal and a length separately is
what makes the geometric conservation law

```math
\sum_{e\in\partial i}\ell_e\,\boldsymbol m_e = \boldsymbol 0
```

hold *exactly* on a warped panel: each side contributes the in-plane edge vector
rotated by a quarter turn, and those sum to the rotation of a closed polygon,
which is zero regardless of how badly the quadrilateral is warped. Deriving the
sense from the panel winding rather than from a per-edge outward test is
essential — the latter misfires on skewed panels.

Each panel carries its own normal for a shared edge, as in Mughal's scheme;
conservation across the edge is then maintained by rotating the donor flux with
`edge_rotation`, which maps components from the second panel's basis into the
first panel's.
"""
struct SurfaceMetrics{T}
    tangent::Matrix{T}
    binormal::Matrix{T}
    edge_length::Vector{T}
    edge_midpoint::Matrix{T}
    edge_area_vector::Array{T, 3}
    edge_rotation::Array{T, 3}
    cell_area::Vector{T}
end

raw"""
    build_surface_metrics(mesh, topology)

Compute the [`SurfaceMetrics`](@ref) of a mesh.

The local basis is built from the panel's own first edge, ``\hat e_1 =
\widehat{P_n(v_2-v_1)}``, rather than by projecting a global axis. A projected
global axis degenerates exactly where the surface normal aligns with it — at the
stem and the transom of a hull — and the usual remedy of picking the
least-aligned axis introduces a discrete branch that would break
differentiability.
"""
function build_surface_metrics(mesh::Mesh, topology::SurfaceTopology)
    element_type = promote_type(eltype(mesh.vertices), eltype(mesh.centers),
        eltype(mesh.normals), eltype(mesh.areas))
    tangent = zeros(element_type, mesh.nfaces, 3)
    binormal = zeros(element_type, mesh.nfaces, 3)
    for panel in 1:mesh.nfaces
        corners = mesh.faces[panel, :] .+ 1
        edge_vector = @view(mesh.vertices[corners[2], :]) .-
                      @view(mesh.vertices[corners[1], :])
        normal = @view mesh.normals[panel, :]
        projected = edge_vector .- dot(edge_vector, normal) .* normal
        length_of = norm(projected)
        axis_one = length_of > 0 ? projected ./ length_of : _fallback_tangent(normal)
        tangent[panel, :] .= axis_one
        binormal[panel, :] .= cross(normal, axis_one)
    end

    # Winding sense of each panel in its own basis: the sign of the projected
    # polygon's signed area. This fixes, once per panel, whether a quarter turn
    # of an edge vector points out of the panel or into it.
    winding = ones(element_type, mesh.nfaces)
    for panel in 1:mesh.nfaces
        corners = mesh.faces[panel, :] .+ 1
        local_x = ntuple(k -> dot(@view(mesh.vertices[corners[k], :]),
                @view(tangent[panel, :])), 4)
        local_y = ntuple(k -> dot(@view(mesh.vertices[corners[k], :]),
                @view(binormal[panel, :])), 4)
        signed_area = zero(element_type)
        for k in 1:4
            next = mod1(k + 1, 4)
            signed_area += local_x[k] * local_y[next] - local_x[next] * local_y[k]
        end
        winding[panel] = signed_area >= 0 ? one(element_type) : -one(element_type)
    end

    edge_length = zeros(element_type, topology.nedges)
    edge_midpoint = zeros(element_type, topology.nedges, 3)
    edge_area_vector = zeros(element_type, topology.nedges, 2, 2)
    edge_rotation = zeros(element_type, topology.nedges, 2, 2)

    for edge in 1:topology.nedges
        for slot in 1:2
            panel = topology.edge_cells[edge, slot]
            panel == 0 && continue
            side = topology.edge_local_side[edge, slot]
            # Index the two corners directly rather than materialising the whole
            # face row: `mesh.faces[panel, :] .+ 1` allocates a vector for every
            # edge slot, and this loop runs twice per edge.
            start = row3(mesh.vertices, mesh.faces[panel, side] + 1)
            finish = row3(mesh.vertices, mesh.faces[panel, mod1(side + 1, 4)] + 1)
            along = finish - start
            panel_tangent = row3(tangent, panel)
            panel_binormal = row3(binormal, panel)
            # In-plane edge vector, turned a quarter turn. Summed over the
            # panel's sides this telescopes to zero exactly.
            along_x = dot(along, panel_tangent)
            along_y = dot(along, panel_binormal)
            edge_area_vector[edge, slot, 1] = winding[panel] * along_y
            edge_area_vector[edge, slot, 2] = -winding[panel] * along_x
            if slot == 1
                edge_length[edge] = norm(along)
                edge_midpoint[edge, :] .= (start .+ finish) ./ 2
            end
        end

        left = topology.edge_cells[edge, 1]
        right = topology.edge_cells[edge, 2]
        if right == 0
            edge_rotation[edge, 1, 1] = one(element_type)
            edge_rotation[edge, 2, 2] = one(element_type)
        else
            left_tangent = @view tangent[left, :]
            left_binormal = @view binormal[left, :]
            right_tangent = @view tangent[right, :]
            right_binormal = @view binormal[right, :]
            # Components of the right panel's basis resolved in the left panel's.
            edge_rotation[edge, 1, 1] = dot(right_tangent, left_tangent)
            edge_rotation[edge, 1, 2] = dot(right_binormal, left_tangent)
            edge_rotation[edge, 2, 1] = dot(right_tangent, left_binormal)
            edge_rotation[edge, 2, 2] = dot(right_binormal, left_binormal)
        end
    end

    return SurfaceMetrics(tangent, binormal, edge_length, edge_midpoint,
        edge_area_vector, edge_rotation, collect(element_type, mesh.areas))
end

# Only reached for a degenerate panel, where the first edge is parallel to the
# normal. Any tangent will do; pick one deterministically.
function _fallback_tangent(normal)
    reference = abs(normal[1]) < 0.9 ? (1.0, 0.0, 0.0) : (0.0, 1.0, 0.0)
    candidate = cross(collect(reference), normal)
    return candidate ./ norm(candidate)
end

"""
    edge_area_vector_for(topology, metrics, edge, panel)

Outward edge normal times edge length for `edge` as seen by `panel`, in that
panel's own basis. This is the quantity a finite-volume flux is weighted by.
"""
function edge_area_vector_for(topology::SurfaceTopology, metrics::SurfaceMetrics,
        edge::Integer, panel::Integer)
    slot = _edge_slot(topology, edge, panel)
    return (metrics.edge_area_vector[edge, slot, 1],
        metrics.edge_area_vector[edge, slot, 2])
end

"""
    edge_normal_for(topology, metrics, edge, panel)

Outward in-surface unit normal of `edge` as seen by `panel`, in that panel's own
basis.
"""
function edge_normal_for(topology::SurfaceTopology, metrics::SurfaceMetrics,
        edge::Integer, panel::Integer)
    area_vector = edge_area_vector_for(topology, metrics, edge, panel)
    magnitude = hypot(area_vector[1], area_vector[2])
    magnitude > 0 || return (one(magnitude), zero(magnitude))
    return (area_vector[1] / magnitude, area_vector[2] / magnitude)
end

function _edge_slot(topology::SurfaceTopology, edge::Integer, panel::Integer)
    topology.edge_cells[edge, 1] == panel && return 1
    topology.edge_cells[edge, 2] == panel && return 2
    throw(ArgumentError("panel $(panel) does not own edge $(edge)"))
end
