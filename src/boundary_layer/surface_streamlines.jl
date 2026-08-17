using LinearAlgebra: dot, norm

raw"""
    SurfaceStreamline

A streamline traced over a panelized surface. `points` is an ``N\times3``
matrix ordered in the requested integration direction and `speed` contains the
corresponding tangential speed at every point.
"""
struct SurfaceStreamline{P,S}
    points::P
    speed::S
end

function _nearest_panel_indices(
    mesh::Mesh,
    point::AbstractVector,
    neighbor_count::Integer,
)
    distances_squared = Vector{promote_type(
        eltype(mesh.centers),
        eltype(point),
    )}(undef, mesh.nfaces)
    for panel in 1:mesh.nfaces
        dx = mesh.centers[panel, 1] - point[1]
        dy = mesh.centers[panel, 2] - point[2]
        dz = mesh.centers[panel, 3] - point[3]
        distances_squared[panel] = dx * dx + dy * dy + dz * dz
    end
    count = min(neighbor_count, mesh.nfaces)
    return partialsortperm(distances_squared, 1:count), distances_squared
end

function _interpolate_surface_velocity(
    mesh::Mesh,
    velocity::AbstractMatrix,
    point::AbstractVector,
    neighbor_count::Integer,
)
    panels, distances_squared = _nearest_panel_indices(
        mesh,
        point,
        neighbor_count,
    )
    T = promote_type(
        eltype(mesh.centers),
        eltype(mesh.normals),
        eltype(velocity),
        eltype(point),
    )
    local_length_squared = sum(mesh.areas[panels]) / length(panels)
    regularization = max(local_length_squared / 16, eps(Float64))
    center = zeros(T, 3)
    normal = zeros(T, 3)
    local_velocity = zeros(T, 3)
    weight_sum = zero(T)
    for panel in panels
        weight = inv(distances_squared[panel] + regularization)
        weight_sum += weight
        for axis in 1:3
            center[axis] += weight * mesh.centers[panel, axis]
            normal[axis] += weight * mesh.normals[panel, axis]
            local_velocity[axis] += weight * velocity[panel, axis]
        end
    end
    center ./= weight_sum
    local_velocity ./= weight_sum
    normal_norm = norm(normal)
    normal_norm > 0 || throw(ArgumentError(
        "neighboring panel normals cancel at the interpolation point",
    ))
    normal ./= normal_norm
    surface_point = collect(point) .- dot(collect(point) .- center, normal) .* normal
    tangential_velocity = local_velocity .- dot(local_velocity, normal) .* normal
    distance_to_surface = sqrt(minimum(distances_squared))
    return surface_point, tangential_velocity, distance_to_surface
end

raw"""
    surface_stagnation_panel(mesh, velocity; panel_mask=nothing)

Return the lowest-speed panel in a supplied search region. The velocity must
contain the complete body-relative surface velocity, including the uniform
incoming stream. The result contains `panel`, `point`, `velocity`, and `speed`.
"""
function surface_stagnation_panel(
    mesh::Mesh,
    velocity::AbstractMatrix;
    panel_mask::Union{Nothing,AbstractVector{Bool}}=nothing,
)
    size(velocity) == (mesh.nfaces, 3) || throw(DimensionMismatch(
        "velocity must have size (mesh.nfaces, 3)",
    ))
    mask = isnothing(panel_mask) ? trues(mesh.nfaces) : BitVector(panel_mask)
    length(mask) == mesh.nfaces || throw(DimensionMismatch(
        "panel_mask must contain one value per mesh face",
    ))
    panels = findall(mask)
    isempty(panels) && throw(ArgumentError(
        "panel_mask must select at least one panel",
    ))
    speeds = [norm(@view velocity[panel, :]) for panel in panels]
    panel = panels[argmin(speeds)]
    return (
        panel,
        point=collect(@view mesh.centers[panel, :]),
        velocity=collect(@view velocity[panel, :]),
        speed=norm(@view velocity[panel, :]),
    )
end

raw"""
    trace_surface_streamline(
        mesh,
        velocity,
        seed;
        step_size=nothing,
        max_steps=500,
        neighbor_count=8,
        direction=1,
        minimum_speed=nothing,
        maximum_surface_distance=nothing,
    )

Trace a streamline over a quadrilateral panel surface using inverse-distance
interpolation of neighboring panel-center velocities. Every integration step
is projected back to the local tangent surface. `seed` may be a panel index or
a three-component point.

The input must be the complete body-relative tangential velocity. For a ship
whose positive ``x`` axis points toward the bow, forward motion therefore has
an incoming-stream contribution ``-U\boldsymbol{e}_x``. Plotting only the
perturbation-potential gradient does not reveal the bow stagnation topology.
"""
function trace_surface_streamline(
    mesh::Mesh,
    velocity::AbstractMatrix,
    seed;
    step_size=nothing,
    max_steps::Integer=500,
    neighbor_count::Integer=8,
    direction::Integer=1,
    minimum_speed=nothing,
    maximum_surface_distance=nothing,
)
    size(velocity) == (mesh.nfaces, 3) || throw(DimensionMismatch(
        "velocity must have size (mesh.nfaces, 3)",
    ))
    max_steps > 0 || throw(ArgumentError("max_steps must be positive"))
    neighbor_count > 0 || throw(ArgumentError(
        "neighbor_count must be positive",
    ))
    direction in (-1, 1) || throw(ArgumentError(
        "direction must be -1 or 1",
    ))
    point = if seed isa Integer
        1 <= seed <= mesh.nfaces || throw(BoundsError(mesh.areas, seed))
        collect(@view mesh.centers[seed, :])
    else
        length(seed) == 3 || throw(DimensionMismatch(
            "seed must contain three coordinates",
        ))
        collect(seed)
    end
    characteristic_length = sqrt(sum(mesh.areas) / mesh.nfaces)
    ds = isnothing(step_size) ? characteristic_length / 3 : step_size
    ds > 0 || throw(ArgumentError("step_size must be positive"))
    speed_floor = isnothing(minimum_speed) ?
        sqrt(eps(Float64)) * maximum(norm(@view velocity[i, :]) for i in 1:mesh.nfaces) :
        minimum_speed
    speed_floor >= 0 || throw(ArgumentError(
        "minimum_speed must be nonnegative",
    ))
    distance_limit = isnothing(maximum_surface_distance) ?
        5 * characteristic_length : maximum_surface_distance
    distance_limit > 0 || throw(ArgumentError(
        "maximum_surface_distance must be positive",
    ))

    T = promote_type(eltype(point), eltype(velocity), eltype(mesh.centers))
    lower_bound = [minimum(mesh.vertices[:, axis]) for axis in 1:3]
    upper_bound = [maximum(mesh.vertices[:, axis]) for axis in 1:3]
    bound_padding = 2 * ds
    point_rows = NTuple{3,T}[]
    speeds = T[]
    surface_point, local_velocity, distance = _interpolate_surface_velocity(
        mesh,
        velocity,
        point,
        neighbor_count,
    )
    distance <= distance_limit || throw(ArgumentError(
        "seed is farther than maximum_surface_distance from the panel surface",
    ))
    for _ in 1:max_steps
        local_speed = norm(local_velocity)
        push!(point_rows, (surface_point[1], surface_point[2], surface_point[3]))
        push!(speeds, local_speed)
        local_speed > speed_floor || break
        candidate = surface_point .+
            direction * ds .* local_velocity ./ local_speed
        next_point, next_velocity, next_distance =
            _interpolate_surface_velocity(
                mesh,
                velocity,
                candidate,
                neighbor_count,
            )
        next_distance <= distance_limit || break
        any(next_point .< lower_bound .- bound_padding) && break
        any(next_point .> upper_bound .+ bound_padding) && break
        norm(next_point .- surface_point) > eps(Float64) * characteristic_length ||
            break
        surface_point = next_point
        local_velocity = next_velocity
    end
    points = Matrix{T}(undef, length(point_rows), 3)
    for (row, values) in enumerate(point_rows), axis in 1:3
        points[row, axis] = values[axis]
    end
    return SurfaceStreamline(points, speeds)
end

raw"""
    trace_surface_streamlines(mesh, velocity, seeds; kwargs...)

Trace one surface streamline from each entry of `seeds`. Keyword arguments are
forwarded to [`trace_surface_streamline`](@ref).
"""
function trace_surface_streamlines(
    mesh::Mesh,
    velocity::AbstractMatrix,
    seeds;
    kwargs...,
)
    return [
        trace_surface_streamline(mesh, velocity, seed; kwargs...)
        for seed in seeds
    ]
end
