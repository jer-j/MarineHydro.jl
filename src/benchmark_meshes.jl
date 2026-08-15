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
    quadrilaterals = Matrix{Float64}[]
    for j in 1:(size(patch, 2) - 1), i in 1:(size(patch, 1) - 1)
        push!(quadrilaterals, Float64[
            patch[i, j, 1] patch[i, j, 2] patch[i, j, 3]
            patch[i + 1, j, 1] patch[i + 1, j, 2] patch[i + 1, j, 3]
            patch[i + 1, j + 1, 1] patch[i + 1, j + 1, 2] patch[i + 1, j + 1, 3]
            patch[i, j + 1, 1] patch[i, j + 1, 2] patch[i, j + 1, 3]
        ])
    end
    return quadrilaterals
end

function _orient_patch!(quadrilaterals, interior_point)
    orientation_score = 0.0
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
