using LinearAlgebra: cross, dot, norm

raw"""
    StructuredPanelGrid

A quadrilateral [`Mesh`](@ref) together with its longitudinal surface-strip
topology. Each entry of `strips` contains panel indices ordered from bow to
stern. The topology is retained for surface-marching methods such as integral
boundary-layer solvers, while `mesh` remains directly usable by the boundary
element routines.
"""
struct StructuredPanelGrid{M,S}
    mesh::M
    strips::S
end

raw"""
    read_tecplot_structured_surfaces(path)

Read one or more Tecplot ASCII `ZONE I=... J=... F=POINT` surface patches.
Each returned patch is an ``I \times J \times 3`` `Array{Float64,3}` whose last
index contains ``(x,y,z)``.

This deliberately small reader targets the structured point files distributed
for the Gothenburg ship-hydrodynamics workshops. It rejects incomplete zones
and does not interpret cell-centered or block-formatted Tecplot data.
"""
function read_tecplot_structured_surfaces(path::AbstractString)
    lines = readlines(path)
    patches = Array{Float64,3}[]
    line_index = 1

    while line_index <= length(lines)
        line = strip(lines[line_index])
        if !startswith(lowercase(line), "zone")
            line_index += 1
            continue
        end

        i_match = match(r"(?i)\bi\s*=\s*(\d+)", line)
        j_match = match(r"(?i)\bj\s*=\s*(\d+)", line)
        isnothing(i_match) && throw(ArgumentError(
            "Tecplot zone on line $(line_index) has no I dimension",
        ))
        isnothing(j_match) && throw(ArgumentError(
            "Tecplot zone on line $(line_index) has no J dimension",
        ))
        ni = parse(Int, something(i_match.captures[1]))
        nj = parse(Int, something(j_match.captures[1]))
        ni >= 2 || throw(ArgumentError("Tecplot I dimension must be at least 2"))
        nj >= 2 || throw(ArgumentError("Tecplot J dimension must be at least 2"))

        expected_values = 3 * ni * nj
        values = Float64[]
        line_index += 1
        while length(values) < expected_values && line_index <= length(lines)
            data_line = strip(lines[line_index])
            startswith(lowercase(data_line), "zone") && throw(ArgumentError(
                "Tecplot zone is incomplete before line $(line_index)",
            ))
            normalized_line = replace(
                data_line,
                'D' => 'E',
                'd' => 'e',
                ',' => ' ',
            )
            append!(values, parse.(Float64, split(normalized_line)))
            line_index += 1
        end
        length(values) == expected_values || throw(ArgumentError(
            "Tecplot zone contains $(length(values)) values; expected $(expected_values)",
        ))

        coordinate_first = reshape(values, 3, ni, nj)
        push!(patches, copy(permutedims(coordinate_first, (2, 3, 1))))
    end

    isempty(patches) && throw(ArgumentError(
        "no structured Tecplot zones were found in $(path)",
    ))
    return patches
end

raw"""
    read_plot3d_structured_surfaces(path)

Read a formatted, block-ordered Plot3D surface file. The first value is the
number of blocks, followed by one ``(n_i,n_j,n_k)`` triplet per block and the
block-ordered ``x``, ``y``, and ``z`` coordinates. Surface blocks must have
``n_k=1``. Each returned patch has size ``n_i\times n_j\times3``.

This reader covers the public DTMB 5415 ``5415_static.net`` surface supplied
by NMRI for the Gothenburg 2000 workshop. Fortran ``D`` exponents are accepted.
"""
function read_plot3d_structured_surfaces(path::AbstractString)
    normalized = replace(read(path, String), 'D' => 'E', 'd' => 'e')
    tokens = split(normalized)
    isempty(tokens) && throw(ArgumentError("Plot3D file is empty"))
    block_count = parse(Int, tokens[1])
    block_count >= 1 || throw(ArgumentError(
        "Plot3D block count must be positive",
    ))
    dimension_end = 1 + 3block_count
    length(tokens) >= dimension_end || throw(ArgumentError(
        "Plot3D file ends before the block dimensions are complete",
    ))

    dimensions = ntuple(block_count) do block
        offset = 1 + 3(block - 1)
        dimensions_for_block = (
            parse(Int, tokens[offset + 1]),
            parse(Int, tokens[offset + 2]),
            parse(Int, tokens[offset + 3]),
        )
        dimensions_for_block[1] >= 2 || throw(ArgumentError(
            "Plot3D surface block $(block) must have at least two i points",
        ))
        dimensions_for_block[2] >= 2 || throw(ArgumentError(
            "Plot3D surface block $(block) must have at least two j points",
        ))
        dimensions_for_block[3] == 1 || throw(ArgumentError(
            "Plot3D block $(block) is volumetric; expected n_k=1",
        ))
        dimensions_for_block
    end

    token_index = dimension_end + 1
    patches = Array{Float64,3}[]
    for (block, (ni, nj, _)) in enumerate(dimensions)
        point_count = ni * nj
        coordinate_count = 3point_count
        token_index + coordinate_count - 1 <= length(tokens) || throw(ArgumentError(
            "Plot3D block $(block) is incomplete",
        ))
        values = parse.(
            Float64,
            @view(tokens[token_index:(token_index + coordinate_count - 1)]),
        )
        coordinate_first = reshape(values, point_count, 3)
        patch = Array{Float64}(undef, ni, nj, 3)
        for coordinate in 1:3
            patch[:, :, coordinate] .= reshape(
                @view(coordinate_first[:, coordinate]),
                ni,
                nj,
            )
        end
        push!(patches, patch)
        token_index += coordinate_count
    end
    token_index == length(tokens) + 1 || throw(ArgumentError(
        "Plot3D file contains $(length(tokens) - token_index + 1) trailing values",
    ))
    return patches
end

function _endpoint_indices(point_count::Integer, requested_count::Integer)
    point_count >= 2 || throw(ArgumentError("a structured direction needs at least 2 points"))
    requested_count >= 2 || throw(ArgumentError("target_shape entries must be at least 2"))
    count = min(point_count, requested_count)
    return unique(round.(Int, range(1, point_count; length=count)))
end

function _transform_gothenburg_patch(
    patch::AbstractArray{<:Real,3},
    length_scale::Real,
    target_shape::Union{Nothing,Tuple{Int,Int}},
)
    size(patch, 3) == 3 || throw(DimensionMismatch(
        "structured patch must have three coordinate components",
    ))
    if isnothing(target_shape)
        i_indices = collect(axes(patch, 1))
        j_indices = collect(axes(patch, 2))
    else
        i_indices = _endpoint_indices(size(patch, 1), target_shape[1])
        j_indices = _endpoint_indices(size(patch, 2), target_shape[2])
    end

    transformed = Array{Float64}(undef, length(i_indices), length(j_indices), 3)
    for (new_j, old_j) in enumerate(j_indices), (new_i, old_i) in enumerate(i_indices)
        # The workshop coordinates have x positive aft and y positive to
        # starboard. MarineHydro's maneuvering convention has x toward the bow
        # and y toward port.
        transformed[new_i, new_j, 1] = -length_scale * patch[old_i, old_j, 1]
        transformed[new_i, new_j, 2] = -length_scale * patch[old_i, old_j, 2]
        transformed[new_i, new_j, 3] = length_scale * patch[old_i, old_j, 3]
    end
    return transformed
end

raw"""
    structured_sectional_area_curve(patch; mirror=true)

Calculate the transverse sectional-area curve of a structured half-hull patch.
The first patch index must run longitudinally and the second must trace each
section in the transverse-vertical plane. The surface trace is closed through
the coordinate origin. This is appropriate for the Gothenburg single-block
files, whose section endpoints lie on the centerplane and design waterline.

For a mirrored half section,

```math
A_i=2\left|\frac{1}{2}\sum_j
\left(y_{i,j}z_{i,j+1}-y_{i,j+1}z_{i,j}\right)\right|.
```

Set `mirror=false` to return the half-sectional area. The result contains
vectors named `x` and `area`.
"""
function structured_sectional_area_curve(
    patch::AbstractArray{<:Real,3};
    mirror::Bool=true,
)
    size(patch, 1) >= 2 || throw(DimensionMismatch(
        "structured patch must contain at least two longitudinal stations",
    ))
    size(patch, 2) >= 2 || throw(DimensionMismatch(
        "each structured section must contain at least two points",
    ))
    size(patch, 3) == 3 || throw(DimensionMismatch(
        "structured patch must have three coordinate components",
    ))

    T = promote_type(Float64, eltype(patch))
    station_x = Vector{T}(undef, size(patch, 1))
    section_area = Vector{T}(undef, size(patch, 1))
    area_multiplier = mirror ? one(T) : one(T) / 2

    for (output_index, station) in enumerate(axes(patch, 1))
        x_sum = zero(T)
        twice_half_area = zero(T)
        for point in axes(patch, 2)
            x_sum += patch[station, point, 1]
            if point < lastindex(patch, 2)
                y_current = patch[station, point, 2]
                z_current = patch[station, point, 3]
                y_next = patch[station, point + 1, 2]
                z_next = patch[station, point + 1, 3]
                twice_half_area += y_current * z_next - y_next * z_current
            end
        end
        station_x[output_index] = x_sum / size(patch, 2)
        section_area[output_index] = area_multiplier * abs(twice_half_area)
    end
    return (x=station_x, area=section_area)
end

raw"""
    gothenburg_maximum_section(
        paths;
        length_scale=1,
        relative_area_tolerance=1e-3,
    )

Find a Schmitz stern-truncation section from Gothenburg single-block half-hull
files. Coordinates are transformed to MarineHydro's ``x``-forward convention.
The selected section is the aftmost station satisfying

```math
A(x)\geq (1-\epsilon)A_{max},
```

where ``\epsilon`` is `relative_area_tolerance`. This places the cutoff at the
aft end of a numerically near-constant parallel middle body. Set the tolerance
to zero to select the discrete maximum exactly.

The return value contains `x`, `area`, `maximum_x`, `maximum_area`,
`relative_area_tolerance`, and `candidate_count`.
"""
function gothenburg_maximum_section(
    paths::AbstractVector{<:AbstractString};
    length_scale::Real=1,
    relative_area_tolerance::Real=1e-3,
)
    isempty(paths) && throw(ArgumentError("at least one surface file is required"))
    length_scale > 0 || throw(ArgumentError("length_scale must be positive"))
    0 <= relative_area_tolerance < 1 || throw(ArgumentError(
        "relative_area_tolerance must be in [0, 1)",
    ))

    station_x = Float64[]
    section_area = Float64[]
    for path in paths
        for patch in read_tecplot_structured_surfaces(path)
            transformed = _transform_gothenburg_patch(
                patch,
                length_scale,
                nothing,
            )
            curve = structured_sectional_area_curve(transformed)
            append!(station_x, curve.x)
            append!(section_area, curve.area)
        end
    end

    maximum_index = argmax(section_area)
    maximum_area = section_area[maximum_index]
    threshold = (1 - relative_area_tolerance) * maximum_area
    candidate_indices = findall(area -> area >= threshold, section_area)
    selected_index = argmin(index -> station_x[index], candidate_indices)

    return (
        x=station_x[selected_index],
        area=section_area[selected_index],
        maximum_x=station_x[maximum_index],
        maximum_area,
        relative_area_tolerance=Float64(relative_area_tolerance),
        candidate_count=length(candidate_indices),
    )
end

gothenburg_maximum_section(paths::AbstractString...; kwargs...) =
    gothenburg_maximum_section(collect(paths); kwargs...)

function _quad_geometry(vertices::AbstractMatrix{<:Real})
    size(vertices) == (4, 3) || throw(DimensionMismatch(
        "quadrilateral vertices must have size (4, 3)",
    ))
    p1 = @view vertices[1, :]
    p2 = @view vertices[2, :]
    p3 = @view vertices[3, :]
    p4 = @view vertices[4, :]
    first_area_vector = cross(p2 - p1, p3 - p1) / 2
    second_area_vector = cross(p3 - p1, p4 - p1) / 2
    first_area = norm(first_area_vector)
    second_area = norm(second_area_vector)
    area = first_area + second_area
    area_vector = first_area_vector + second_area_vector
    normal_magnitude = norm(area_vector)

    if iszero(area) || iszero(normal_magnitude)
        return nothing
    end

    first_centroid = (p1 + p2 + p3) / 3
    second_centroid = (p1 + p3 + p4) / 3
    center = (
        first_area * first_centroid + second_area * second_centroid
    ) / area
    panel_normal = area_vector / normal_magnitude
    radius = maximum(norm(@view(vertices[index, :]) - center) for index in 1:4)
    return (
        center=collect(center),
        normal=collect(panel_normal),
        area=area,
        radius=radius,
        area_vector=collect(area_vector),
    )
end

function _patch_quadrilaterals(patch::AbstractArray{<:Real,3})
    T = eltype(patch)
    quadrilaterals = Matrix{T}[]
    for j in 1:(size(patch, 2) - 1), i in 1:(size(patch, 1) - 1)
        push!(quadrilaterals, T[
            patch[i, j, 1] patch[i, j, 2] patch[i, j, 3]
            patch[i + 1, j, 1] patch[i + 1, j, 2] patch[i + 1, j, 3]
            patch[i + 1, j + 1, 1] patch[i + 1, j + 1, 2] patch[i + 1, j + 1, 3]
            patch[i, j + 1, 1] patch[i, j + 1, 2] patch[i, j + 1, 3]
        ])
    end
    return quadrilaterals
end

function _orient_patch!(quadrilaterals, interior_point)
    isempty(quadrilaterals) && return quadrilaterals
    T = promote_type(eltype(first(quadrilaterals)), eltype(interior_point))
    orientation_score = zero(T)
    for vertices in quadrilaterals
        geometry = _quad_geometry(vertices)
        isnothing(geometry) && continue
        orientation_score += dot(
            geometry.area_vector,
            geometry.center - interior_point,
        )
    end
    if orientation_score < 0
        for index in eachindex(quadrilaterals)
            quadrilaterals[index] = quadrilaterals[index][[1, 4, 3, 2], :]
        end
    end
    return quadrilaterals
end

function _matrix_from_rows(rows::Vector{NTuple{N,T}}) where {N,T}
    matrix = Matrix{T}(undef, length(rows), N)
    for (row_index, row) in enumerate(rows), column_index in 1:N
        matrix[row_index, column_index] = row[column_index]
    end
    return matrix
end

function _clip_polygon_below_waterline(vertices, waterline)
    output = Vector{Vector{Float64}}()
    isempty(vertices) && return output
    previous = Float64.(vertices[end])
    previous_inside = previous[3] <= waterline
    for current_input in vertices
        current = Float64.(current_input)
        current_inside = current[3] <= waterline
        if current_inside != previous_inside
            fraction = (waterline - previous[3]) / (current[3] - previous[3])
            intersection = previous .+ fraction .* (current .- previous)
            intersection[3] = waterline
            push!(output, intersection)
        end
        current_inside && push!(output, current)
        previous = current
        previous_inside = current_inside
    end
    return output
end

function _quadrilaterals_from_polygon(polygon)
    length(polygon) < 3 && return Matrix{Float64}[]
    if length(polygon) == 3
        return [reduce(vcat, permutedims.(polygon[[1, 2, 3, 3]]))]
    elseif length(polygon) == 4
        return [reduce(vcat, permutedims.(polygon))]
    end
    return [
        reduce(vcat, permutedims.(polygon[[1, index, index + 1, index + 1]]))
        for index in 2:(length(polygon) - 1)
    ]
end

function _mesh_from_quadrilaterals(
    quadrilaterals::AbstractVector;
    interior_point,
    area_tolerance,
)
    isempty(quadrilaterals) && throw(ArgumentError(
        "no immersed nondegenerate panels were generated",
    ))
    _orient_patch!(quadrilaterals, interior_point)
    T = eltype(first(quadrilaterals))
    vertex_rows = NTuple{3,T}[]
    face_rows = NTuple{4,Int}[]
    center_rows = NTuple{3,T}[]
    normal_rows = NTuple{3,T}[]
    areas = T[]
    radii = T[]
    for vertices in quadrilaterals
        geometry = _quad_geometry(vertices)
        if isnothing(geometry) || geometry.area <= area_tolerance
            continue
        end
        first_vertex = length(vertex_rows)
        append!(vertex_rows, [Tuple(row) for row in eachrow(vertices)])
        push!(face_rows, (
            first_vertex,
            first_vertex + 1,
            first_vertex + 2,
            first_vertex + 3,
        ))
        push!(center_rows, Tuple(geometry.center))
        push!(normal_rows, Tuple(geometry.normal))
        push!(areas, geometry.area)
        push!(radii, geometry.radius)
    end
    isempty(face_rows) && throw(ArgumentError(
        "all panels were degenerate at the requested area tolerance",
    ))
    return Mesh(
        _matrix_from_rows(vertex_rows),
        _matrix_from_rows(face_rows),
        _matrix_from_rows(center_rows),
        _matrix_from_rows(normal_rows),
        areas,
        radii,
        length(vertex_rows),
        length(face_rows),
    )
end

raw"""
    structured_surface_waterline_mesh(
        patches;
        coordinate_transform=identity,
        waterline=0,
        target_shape=nothing,
        mirror=false,
        interior_point=nothing,
        area_tolerance=nothing,
    )

Clip structured surface patches to ``z\leq z_w`` and construct an immersed
MarineHydro mesh. `coordinate_transform` maps each source point into
MarineHydro coordinates before clipping. Cells crossing the waterline are
intersected exactly and split into triangular panels when necessary. A
triangle is represented by a four-node panel whose final edge has zero length,
which is supported by the Rankine panel integral.

Set `mirror=true` for a half hull. The input side is reflected through
``y=0`` after clipping. `target_shape=(n_i,n_j)` applies endpoint-preserving
structured downsampling before the coordinate transform.
"""
function structured_surface_waterline_mesh(
    patches::AbstractVector{<:AbstractArray{<:Real,3}};
    coordinate_transform=identity,
    waterline::Real=0,
    target_shape::Union{Nothing,Tuple{Int,Int}}=nothing,
    mirror::Bool=false,
    interior_point=nothing,
    area_tolerance=nothing,
)
    isempty(patches) && throw(ArgumentError("at least one patch is required"))
    quadrilaterals = Matrix{Float64}[]
    transformed_points = Vector{Float64}[]
    for source_patch in patches
        size(source_patch, 3) == 3 || throw(DimensionMismatch(
            "structured patch must have three coordinate components",
        ))
        i_indices = isnothing(target_shape) ? collect(axes(source_patch, 1)) :
            _endpoint_indices(size(source_patch, 1), target_shape[1])
        j_indices = isnothing(target_shape) ? collect(axes(source_patch, 2)) :
            _endpoint_indices(size(source_patch, 2), target_shape[2])
        patch = Array{Float64}(undef, length(i_indices), length(j_indices), 3)
        for (new_j, old_j) in enumerate(j_indices),
            (new_i, old_i) in enumerate(i_indices)
            point = coordinate_transform(collect(@view source_patch[old_i, old_j, :]))
            length(point) == 3 || throw(DimensionMismatch(
                "coordinate_transform must return three coordinates",
            ))
            patch[new_i, new_j, :] .= point
            push!(transformed_points, Float64.(point))
        end
        for j in 1:(size(patch, 2) - 1), i in 1:(size(patch, 1) - 1)
            polygon = _clip_polygon_below_waterline(
                [
                    collect(@view patch[i, j, :]),
                    collect(@view patch[i + 1, j, :]),
                    collect(@view patch[i + 1, j + 1, :]),
                    collect(@view patch[i, j + 1, :]),
                ],
                waterline,
            )
            append!(quadrilaterals, _quadrilaterals_from_polygon(polygon))
        end
    end
    immersed_points = filter(point -> point[3] <= waterline, transformed_points)
    isempty(immersed_points) && throw(ArgumentError(
        "the structured surface does not intersect or lie below the waterline",
    ))
    point_matrix = reduce(vcat, permutedims.(immersed_points))
    spans = [maximum(point_matrix[:, axis]) - minimum(point_matrix[:, axis]) for axis in 1:3]
    interior = if isnothing(interior_point)
        [
            (minimum(point_matrix[:, 1]) + maximum(point_matrix[:, 1])) / 2,
            0.0,
            (minimum(point_matrix[:, 3]) + waterline) / 2,
        ]
    else
        length(interior_point) == 3 || throw(DimensionMismatch(
            "interior_point must contain three coordinates",
        ))
        Float64.(collect(interior_point))
    end
    if mirror
        source_panels = copy(quadrilaterals)
        for vertices in source_panels
            mirrored = copy(vertices[[1, 4, 3, 2], :])
            mirrored[:, 2] .*= -1
            push!(quadrilaterals, mirrored)
        end
    end
    tolerance = isnothing(area_tolerance) ?
        sqrt(eps(Float64)) * maximum(spans)^2 : Float64(area_tolerance)
    tolerance >= 0 || throw(ArgumentError("area_tolerance must be nonnegative"))
    return _mesh_from_quadrilaterals(
        quadrilaterals;
        interior_point=interior,
        area_tolerance=tolerance,
    )
end

structured_surface_waterline_mesh(patch::AbstractArray{<:Real,3}; kwargs...) =
    structured_surface_waterline_mesh([patch]; kwargs...)

raw"""
    read_ascii_stl_waterline_mesh(
        path;
        coordinate_transform=identity,
        waterline=0,
        clustering=(48, 16, 8),
        interior_point=nothing,
        area_tolerance=nothing,
    )

Read an ASCII STL hull, clip its triangles to ``z\leq z_w``, and construct a
waterline-conforming MarineHydro mesh. `coordinate_transform` is applied before
clipping. Vertex clustering on a regular ``(n_x,n_y,n_z)`` lattice provides a
deterministic reduction of very dense public triangulations. Waterline
vertices occupy a separate cluster layer, so the clipped boundary remains at
``z=z_w``.

The input must be plain ASCII STL. Decompress ``.stl.gz`` resources before
calling this function.
"""
function read_ascii_stl_waterline_mesh(
    path::AbstractString;
    coordinate_transform=identity,
    waterline::Real=0,
    clustering::NTuple{3,Int}=(48, 16, 8),
    interior_point=nothing,
    area_tolerance=nothing,
)
    all(count -> count >= 1, clustering) || throw(ArgumentError(
        "clustering entries must be positive",
    ))
    source_vertices = Vector{Float64}[]
    polygons = Vector{Vector{Vector{Float64}}}()
    for line in eachline(path)
        fields = split(strip(line))
        if length(fields) == 4 && lowercase(fields[1]) == "vertex"
            point = parse.(Float64, fields[2:4])
            transformed = coordinate_transform(point)
            length(transformed) == 3 || throw(DimensionMismatch(
                "coordinate_transform must return three coordinates",
            ))
            push!(source_vertices, Float64.(transformed))
            if length(source_vertices) == 3
                clipped = _clip_polygon_below_waterline(
                    source_vertices,
                    waterline,
                )
                length(clipped) >= 3 && push!(polygons, clipped)
                empty!(source_vertices)
            end
        end
    end
    isempty(source_vertices) || throw(ArgumentError(
        "ASCII STL contains an incomplete triangle",
    ))
    isempty(polygons) && throw(ArgumentError(
        "ASCII STL contains no triangles at or below the waterline",
    ))

    all_points = reduce(vcat, [reduce(vcat, permutedims.(polygon)) for polygon in polygons])
    minima = [minimum(all_points[:, axis]) for axis in 1:3]
    maxima = [maximum(all_points[:, axis]) for axis in 1:3]
    spans = maxima .- minima
    cluster_sums = Dict{NTuple{3,Int},Vector{Float64}}()
    cluster_counts = Dict{NTuple{3,Int},Int}()
    polygon_keys = Vector{Vector{NTuple{3,Int}}}(undef, length(polygons))

    function cluster_key(point)
        horizontal = ntuple(2) do axis
            iszero(spans[axis]) && return 0
            clamp(
                floor(Int, clustering[axis] *
                    (point[axis] - minima[axis]) / spans[axis]),
                0,
                clustering[axis] - 1,
            )
        end
        on_waterline = abs(point[3] - waterline) <=
            32eps(Float64) * max(abs(waterline), maximum(spans), 1.0)
        vertical = if on_waterline
            clustering[3]
        elseif iszero(spans[3])
            0
        else
            clamp(
                floor(Int, clustering[3] *
                    (point[3] - minima[3]) / spans[3]),
                0,
                clustering[3] - 1,
            )
        end
        return (horizontal[1], horizontal[2], vertical)
    end

    for (polygon_index, polygon) in enumerate(polygons)
        keys = NTuple{3,Int}[]
        for point in polygon
            key = cluster_key(point)
            push!(keys, key)
            cluster_sums[key] = get(cluster_sums, key, zeros(3)) .+ point
            cluster_counts[key] = get(cluster_counts, key, 0) + 1
        end
        polygon_keys[polygon_index] = keys
    end
    cluster_points = Dict(
        key => value ./ cluster_counts[key]
        for (key, value) in cluster_sums
    )
    quadrilaterals = Matrix{Float64}[]
    seen_faces = Set{Tuple}()
    for keys in polygon_keys
        reduced_keys = unique(keys)
        length(reduced_keys) >= 3 || continue
        signature = Tuple(sort(reduced_keys))
        signature in seen_faces && continue
        push!(seen_faces, signature)
        polygon = [cluster_points[key] for key in reduced_keys]
        append!(quadrilaterals, _quadrilaterals_from_polygon(polygon))
    end
    interior = if isnothing(interior_point)
        [
            (minima[1] + maxima[1]) / 2,
            (minima[2] + maxima[2]) / 2,
            (minima[3] + waterline) / 2,
        ]
    else
        length(interior_point) == 3 || throw(DimensionMismatch(
            "interior_point must contain three coordinates",
        ))
        Float64.(collect(interior_point))
    end
    tolerance = isnothing(area_tolerance) ?
        sqrt(eps(Float64)) * maximum(spans)^2 : Float64(area_tolerance)
    tolerance >= 0 || throw(ArgumentError("area_tolerance must be nonnegative"))
    return _mesh_from_quadrilaterals(
        quadrilaterals;
        interior_point=interior,
        area_tolerance=tolerance,
    )
end

raw"""
    ellipsoid_potential_coefficients(
        semi_axes;
        quadrature_order=96,
    )

Return the three classical potential coefficients of a triaxial ellipsoid
with semi-axes ``(a,b,c)``. The coefficients are

```math
\alpha_i=abc\int_0^\infty
\frac{\mathrm{d}\lambda}
{(a_i^2+\lambda)
\sqrt{(a^2+\lambda)(b^2+\lambda)(c^2+\lambda)}}.
```

They satisfy ``\alpha_x+\alpha_y+\alpha_z=2``. The semi-infinite integral is
mapped to ``t\in(0,1)`` with
``\lambda=[t/(1-t)]^2`` and evaluated by Gauss-Legendre quadrature. The fixed
quadrature rule keeps the calculation differentiable with respect to the
semi-axes.
"""
function ellipsoid_potential_coefficients(
    semi_axes;
    quadrature_order::Integer=96,
)
    length(semi_axes) == 3 || throw(DimensionMismatch(
        "semi_axes must contain exactly three values",
    ))
    all(axis -> axis > 0, semi_axes) || throw(ArgumentError(
        "ellipsoid semi-axes must be positive",
    ))
    quadrature_order >= 8 || throw(ArgumentError(
        "quadrature_order must be at least 8",
    ))

    a, b, c = promote(float.(semi_axes)...)
    axis_squares = (a^2, b^2, c^2)
    axis_product = a * b * c

    # Golub-Welsch nodes and weights on [0, 1]. For Legendre polynomials the
    # first-row eigenvector components squared are the unit-interval weights.
    off_diagonal = [
        index / sqrt(4index^2 - 1)
        for index in 1:(quadrature_order - 1)
    ]
    factorization = eigen(SymTridiagonal(
        zeros(Float64, quadrature_order),
        off_diagonal,
    ))
    nodes = (factorization.values .+ 1) ./ 2
    weights = factorization.vectors[1, :] .^ 2

    coefficients = ntuple(3) do direction
        integral_value = zero(axis_product)
        for quadrature_index in eachindex(nodes)
            t = nodes[quadrature_index]
            one_minus_t = 1 - t
            ratio = t / one_minus_t
            lambda = ratio^2
            jacobian = 2t / one_minus_t^3
            delta = sqrt(prod(axis_square + lambda for axis_square in axis_squares))
            integrand = axis_product * jacobian /
                ((axis_squares[direction] + lambda) * delta)
            integral_value += weights[quadrature_index] * integrand
        end
        integral_value
    end
    return (x=coefficients[1], y=coefficients[2], z=coefficients[3])
end

raw"""
    ellipsoid_added_mass(
        semi_axes;
        density=SETTINGS.rho,
        quadrature_order=96,
    )

Return the translational added masses of a complete triaxial ellipsoid in an
unbounded ideal fluid. If ``\alpha_i`` is returned by
[`ellipsoid_potential_coefficients`](@ref), then

```math
m_i=\rho\frac{4\pi abc}{3}\frac{\alpha_i}{2-\alpha_i}.
```

For a lower half ellipsoid intersected by a rigid plane at ``z=0``, the
same-sign image body completes the ellipsoid and the force integral on the
wetted half is ``m_i/2``.
"""
function ellipsoid_added_mass(
    semi_axes;
    density::Real=SETTINGS.rho,
    quadrature_order::Integer=96,
)
    density > 0 || throw(ArgumentError("density must be positive"))
    length(semi_axes) == 3 || throw(DimensionMismatch(
        "semi_axes must contain exactly three values",
    ))
    a, b, c = promote(float.(semi_axes)...)
    coefficients = ellipsoid_potential_coefficients(
        (a, b, c);
        quadrature_order,
    )
    displaced_mass = density * (4pi / 3) * a * b * c
    added_mass = ntuple(3) do direction
        coefficient = getproperty(coefficients, (:x, :y, :z)[direction])
        displaced_mass * coefficient / (2 - coefficient)
    end
    return (x=added_mass[1], y=added_mass[2], z=added_mass[3])
end

raw"""
    surface_piercing_ellipsoid_panel_grid(
        semi_axes;
        longitudinal_panels=24,
        girth_panels=12,
    )

Generate the wetted lower half of a triaxial ellipsoid centered in the rigid
free-surface plane ``z=0``. The parameterization is

```math
x=a\cos\mu,\qquad
y=b\sin\mu\cos\nu,\qquad
z=-c\sin\mu\sin\nu,
```

with ``0\leq\mu\leq\pi`` from bow to stern and
``0\leq\nu\leq\pi`` around the immersed girth. The returned longitudinal
strips are ordered from positive ``x`` to negative ``x``.
"""
function surface_piercing_ellipsoid_panel_grid(
    semi_axes;
    longitudinal_panels::Integer=24,
    girth_panels::Integer=12,
)
    length(semi_axes) == 3 || throw(DimensionMismatch(
        "semi_axes must contain exactly three values",
    ))
    all(axis -> axis > 0, semi_axes) || throw(ArgumentError(
        "ellipsoid semi-axes must be positive",
    ))
    longitudinal_panels >= 4 || throw(ArgumentError(
        "longitudinal_panels must be at least 4",
    ))
    girth_panels >= 2 || throw(ArgumentError(
        "girth_panels must be at least 2",
    ))

    a, b, c = promote(float.(semi_axes)...)
    T = promote_type(typeof(a), typeof(b), typeof(c), Float64)
    mu = range(0.0, pi; length=longitudinal_panels + 1)
    nu = range(0.0, pi; length=girth_panels + 1)
    patch = Array{T}(undef, length(mu), length(nu), 3)
    for (girth_index, nu_value) in enumerate(nu),
        (longitudinal_index, mu_value) in enumerate(mu)
        sin_mu = sin(mu_value)
        patch[longitudinal_index, girth_index, 1] = a * cos(mu_value)
        patch[longitudinal_index, girth_index, 2] =
            b * sin_mu * cos(nu_value)
        patch[longitudinal_index, girth_index, 3] =
            -c * sin_mu * sin(nu_value)
    end

    quadrilaterals = _patch_quadrilaterals(patch)
    _orient_patch!(quadrilaterals, T[zero(T), zero(T), -c / 2])
    area_tolerance = sqrt(eps(Float64)) * maximum((a, b, c))^2
    vertex_rows = NTuple{3,T}[]
    face_rows = NTuple{4,Int}[]
    center_rows = NTuple{3,T}[]
    normal_rows = NTuple{3,T}[]
    areas = T[]
    radii = T[]
    strips = [Int[] for _ in 1:girth_panels]

    quadrilateral_index = 0
    for girth_index in 1:girth_panels
        for _ in 1:longitudinal_panels
            quadrilateral_index += 1
            vertices = quadrilaterals[quadrilateral_index]
            geometry = _quad_geometry(vertices)
            if isnothing(geometry) || geometry.area <= area_tolerance
                continue
            end
            first_vertex = length(vertex_rows)
            for row in eachrow(vertices)
                push!(vertex_rows, (row[1], row[2], row[3]))
            end
            push!(face_rows, (
                first_vertex,
                first_vertex + 1,
                first_vertex + 2,
                first_vertex + 3,
            ))
            push!(center_rows, Tuple(geometry.center))
            push!(normal_rows, Tuple(geometry.normal))
            push!(areas, geometry.area)
            push!(radii, geometry.radius)
            push!(strips[girth_index], length(face_rows))
        end
    end
    all(!isempty, strips) || error(
        "ellipsoid discretization produced an empty longitudinal strip",
    )
    mesh = Mesh(
        _matrix_from_rows(vertex_rows),
        _matrix_from_rows(face_rows),
        _matrix_from_rows(center_rows),
        _matrix_from_rows(normal_rows),
        areas,
        radii,
        length(vertex_rows),
        length(face_rows),
    )
    return StructuredPanelGrid(mesh, strips)
end

surface_piercing_ellipsoid_mesh(semi_axes; kwargs...) =
    surface_piercing_ellipsoid_panel_grid(semi_axes; kwargs...).mesh

raw"""
    wigley_hull_panel_grid(
        length,
        breadth,
        draft;
        longitudinal_panels=32,
        vertical_panels=8,
    )

Generate the classical parabolic Wigley hull below ``z=0``,

```math
y=\pm\frac{B}{2}\left[1-\left(\frac{2x}{L}\right)^2\right]
\left[1-\left(\frac{z}{T}\right)^2\right],
\qquad -T\leq z\leq0.
```

The two hull sides are discretized independently. Returned surface strips run
from the positive-``x`` bow to the negative-``x`` stern.
"""
function wigley_hull_panel_grid(
    hull_length::Real,
    breadth::Real,
    draft::Real;
    longitudinal_panels::Integer=32,
    vertical_panels::Integer=8,
)
    hull_length > 0 || throw(ArgumentError("length must be positive"))
    breadth > 0 || throw(ArgumentError("breadth must be positive"))
    draft > 0 || throw(ArgumentError("draft must be positive"))
    longitudinal_panels >= 4 || throw(ArgumentError(
        "longitudinal_panels must be at least 4",
    ))
    vertical_panels >= 2 || throw(ArgumentError(
        "vertical_panels must be at least 2",
    ))
    T = promote_type(
        typeof(float(hull_length)),
        typeof(float(breadth)),
        typeof(float(draft)),
    )
    x_values = range(
        hull_length / 2,
        -hull_length / 2;
        length=longitudinal_panels + 1,
    )
    vertical_fraction = range(0.0, 1.0; length=vertical_panels + 1)
    quadrilaterals = Matrix{T}[]
    strip_panel_ranges = UnitRange{Int}[]
    for side in (-1.0, 1.0), vertical_index in 1:vertical_panels
        first_panel = length(quadrilaterals) + 1
        for longitudinal_index in 1:longitudinal_panels
            vertices = Matrix{T}(undef, 4, 3)
            for (corner, (x_index, z_index)) in enumerate((
                (longitudinal_index, vertical_index),
                (longitudinal_index + 1, vertical_index),
                (longitudinal_index + 1, vertical_index + 1),
                (longitudinal_index, vertical_index + 1),
            ))
                x = x_values[x_index]
                eta = vertical_fraction[z_index]
                vertices[corner, 1] = x
                vertices[corner, 2] = side * breadth / 2 *
                    (1 - (2x / hull_length)^2) * (1 - eta^2)
                vertices[corner, 3] = -draft * eta
            end
            if side > 0
                vertices = vertices[[1, 4, 3, 2], :]
            end
            push!(quadrilaterals, vertices)
        end
        push!(strip_panel_ranges, first_panel:length(quadrilaterals))
    end
    mesh = _mesh_from_quadrilaterals(
        quadrilaterals;
        interior_point=T[zero(T), zero(T), -draft / 2],
        area_tolerance=sqrt(eps(Float64)) *
            maximum((hull_length, breadth, draft))^2,
    )
    strips = [collect(panel_range) for panel_range in strip_panel_ranges]
    return StructuredPanelGrid(mesh, strips)
end

wigley_hull_mesh(hull_length, breadth, draft; kwargs...) =
    wigley_hull_panel_grid(hull_length, breadth, draft; kwargs...).mesh

raw"""
    read_gothenburg2010_mesh(
        paths;
        length_scale=1,
        target_shape=nothing,
        mirror=true,
        interior_point=nothing,
        area_tolerance=nothing,
    )

Build a MarineHydro quadrilateral mesh from the structured half-hull surface
files distributed for the Gothenburg workshops.

The source files use nondimensional coordinates with the forward
perpendicular at ``x=-0.5``, aft perpendicular at ``x=0.5``, and positive
``y`` to starboard. This importer applies

```math
(x,y,z)_{MH}=L(-x,-y,z)_{Gothenburg},
```

so MarineHydro receives ``x`` toward the bow and ``y`` toward port.
`length_scale` is the desired ``L_{pp}``. By default all source points are
retained. Set `target_shape=(n_i,n_j)` for endpoint-preserving structured
downsampling. `mirror=true` reflects the supplied starboard half hull to form
the full wetted surface.

Normals are oriented globally on each structured patch using an interior
point. The automatic interior point is the center of the transformed bounding
box at ``y=0``. Pass a known interior point for unusually non-star-shaped
geometries.
"""
function _read_gothenburg2010_panel_grid(
    paths::AbstractVector{<:AbstractString};
    length_scale::Real=1,
    target_shape::Union{Nothing,Tuple{Int,Int}}=nothing,
    mirror::Bool=true,
    interior_point=nothing,
    area_tolerance=nothing,
)
    isempty(paths) && throw(ArgumentError("at least one surface file is required"))
    length_scale > 0 || throw(ArgumentError("length_scale must be positive"))

    patches = Array{Float64,3}[]
    for path in paths
        append!(
            patches,
            _transform_gothenburg_patch.(
                read_tecplot_structured_surfaces(path),
                Ref(length_scale),
                Ref(target_shape),
            ),
        )
    end

    all_points = reduce(vcat, [reshape(patch, :, 3) for patch in patches])
    if isnothing(interior_point)
        interior = [
            (minimum(all_points[:, 1]) + maximum(all_points[:, 1])) / 2,
            0.0,
            (minimum(all_points[:, 3]) + maximum(all_points[:, 3])) / 2,
        ]
    else
        length(interior_point) == 3 || throw(DimensionMismatch(
            "interior_point must contain three coordinates",
        ))
        interior = Float64.(collect(interior_point))
    end

    tolerance = if isnothing(area_tolerance)
        spans = [maximum(all_points[:, axis]) - minimum(all_points[:, axis]) for axis in 1:3]
        sqrt(eps(Float64)) * maximum(spans)^2
    else
        area_tolerance >= 0 || throw(ArgumentError("area_tolerance must be nonnegative"))
        Float64(area_tolerance)
    end

    vertex_rows = NTuple{3,Float64}[]
    face_rows = NTuple{4,Int}[]
    center_rows = NTuple{3,Float64}[]
    normal_rows = NTuple{3,Float64}[]
    areas = Float64[]
    radii = Float64[]

    function append_quadrilateral!(vertices)
        geometry = _quad_geometry(vertices)
        if isnothing(geometry) || geometry.area <= tolerance
            return nothing
        end
        first_vertex = length(vertex_rows)
        for row in eachrow(vertices)
            push!(vertex_rows, (row[1], row[2], row[3]))
        end
        push!(face_rows, (
            first_vertex,
            first_vertex + 1,
            first_vertex + 2,
            first_vertex + 3,
        ))
        push!(center_rows, (
            geometry.center[1],
            geometry.center[2],
            geometry.center[3],
        ))
        push!(normal_rows, (
            geometry.normal[1],
            geometry.normal[2],
            geometry.normal[3],
        ))
        push!(areas, geometry.area)
        push!(radii, geometry.radius)
        return length(face_rows)
    end

    strip_records = NamedTuple[]
    for patch in patches
        quadrilaterals = _patch_quadrilaterals(patch)
        _orient_patch!(quadrilaterals, interior)
        quadrilateral_index = 0
        transverse_cells = size(patch, 2) - 1
        for transverse_index in 1:transverse_cells
            negative_side_indices = Int[]
            positive_side_indices = Int[]
            for _ in 1:(size(patch, 1) - 1)
                quadrilateral_index += 1
                vertices = quadrilaterals[quadrilateral_index]
                negative_index = append_quadrilateral!(vertices)
                isnothing(negative_index) || push!(negative_side_indices, negative_index)
                if mirror
                    mirrored = copy(vertices[[1, 4, 3, 2], :])
                    mirrored[:, 2] .*= -1
                    positive_index = append_quadrilateral!(mirrored)
                    isnothing(positive_index) || push!(positive_side_indices, positive_index)
                end
            end
            isempty(negative_side_indices) || push!(strip_records, (
                side=-1,
                transverse_cells,
                transverse_index,
                panel_indices=negative_side_indices,
            ))
            if mirror && !isempty(positive_side_indices)
                push!(strip_records, (
                    side=1,
                    transverse_cells,
                    transverse_index,
                    panel_indices=positive_side_indices,
                ))
            end
        end
    end

    isempty(face_rows) && throw(ArgumentError(
        "all structured cells were degenerate at the requested tolerance",
    ))
    mesh = Mesh(
        _matrix_from_rows(vertex_rows),
        _matrix_from_rows(face_rows),
        _matrix_from_rows(center_rows),
        _matrix_from_rows(normal_rows),
        areas,
        radii,
        length(vertex_rows),
        length(face_rows),
    )

    # Bow and stern source files commonly carry separate patches with the same
    # transverse topology. Grouping on side and transverse band stitches those
    # patches into one marching strip. Sorting by decreasing x follows the
    # double-body inflow from bow to stern in MarineHydro coordinates.
    strip_keys = unique(
        (record.side, record.transverse_cells, record.transverse_index)
        for record in strip_records
    )
    strips = Vector{Vector{Int}}(undef, length(strip_keys))
    for (strip_index, key) in enumerate(strip_keys)
        indices = reduce(
            vcat,
            (
                record.panel_indices
                for record in strip_records
                if (record.side, record.transverse_cells, record.transverse_index) == key
            ),
        )
        sort!(indices; by=index -> mesh.centers[index, 1], rev=true)
        strips[strip_index] = indices
    end
    return StructuredPanelGrid(mesh, strips)
end

raw"""
    read_gothenburg2010_panel_grid(paths; kwargs...)

Read Gothenburg structured hull patches while preserving longitudinal panel
strips for surface-marching calculations. Geometry transformation, mirroring,
normal orientation, downsampling, and degeneracy handling are identical to
[`read_gothenburg2010_mesh`](@ref).
"""
read_gothenburg2010_panel_grid(
    paths::AbstractVector{<:AbstractString};
    kwargs...,
) = _read_gothenburg2010_panel_grid(paths; kwargs...)

read_gothenburg2010_panel_grid(paths::AbstractString...; kwargs...) =
    read_gothenburg2010_panel_grid(collect(paths); kwargs...)

function read_gothenburg2010_mesh(
    paths::AbstractVector{<:AbstractString};
    kwargs...,
)
    return read_gothenburg2010_panel_grid(paths; kwargs...).mesh
end

read_gothenburg2010_mesh(paths::AbstractString...; kwargs...) =
    read_gothenburg2010_mesh(collect(paths); kwargs...)

raw"""
    mesh_surface_area(mesh)

Return the sum of panel areas.
"""
mesh_surface_area(mesh::Mesh) = sum(mesh.areas)

raw"""
    mesh_signed_volume(mesh; reference=(0, 0, 0))

Estimate enclosed volume by the divergence theorem,

```math
\nabla=\frac{1}{3}\int_S
(\boldsymbol{x}-\boldsymbol{x}_{ref})\mathbin{\boldsymbol{\cdot}}
\boldsymbol{n}\,\mathrm{d}S.
```

Normals must point out of the body. For a surface-piercing hull whose omitted
closure lies in ``z=0``, choose a reference point in that plane. The closure
then contributes exactly zero to the integral.
"""
function mesh_signed_volume(mesh::Mesh; reference=(0, 0, 0))
    length(reference) == 3 || throw(DimensionMismatch(
        "reference must contain three coordinates",
    ))
    offsets = mesh.centers .- reshape(collect(reference), 1, 3)
    radial_normal_component = vec(sum(offsets .* mesh.normals; dims=2))
    return sum(radial_normal_component .* mesh.areas) / 3
end
