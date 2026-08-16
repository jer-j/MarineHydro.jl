using LinearAlgebra: cross, dot, norm
import NonlinearSolve

using NonlinearSolve: NonlinearProblem, SimpleNewtonRaphson

raw"""
    ThreeDimensionalClosure(; kwargs...)

Parameters for the incompressible, fully turbulent three-dimensional integral
boundary layer.

The three-dimensional structure — the momentum-loss tensor, the crossflow
parameter and the skewed wall shear — follows the local-Cartesian formulation
of Mughal, as used by Lokatt and Eller. The scalar closures for the kinetic
energy shape factor, skin friction, slip velocity, equilibrium shear and
dissipation are Drela's incompressible turbulent set, which is the same
lineage and which Lokatt and Eller also adopt.

# Keywords

- `initial_shape_factor`: kinematic shape factor ``H_k`` at the first station.
- `separation_shape_factor`: ``H_k`` at which the attached-flow closure is
  declared invalid. Unlike the two-dimensional case this is a *streamwise*
  criterion; genuine three-dimensional separation occurs where skin-friction
  lines converge and is diagnosed separately from the crossflow angle.
- `minimum_shape_factor`: lower bound enforced on ``H_k`` for closure validity.
- `maximum_crossflow_angle`: bound on ``|\beta_w|`` in radians. The crossflow
  closure is a small-to-moderate-angle correlation and is not meant to be
  extrapolated to fully three-dimensional separated flow.
"""
struct ThreeDimensionalClosure{T}
    initial_shape_factor::T
    separation_shape_factor::T
    minimum_shape_factor::T
    maximum_crossflow_angle::T
end

function ThreeDimensionalClosure(; initial_shape_factor = 1.4,
        separation_shape_factor = 2.4, minimum_shape_factor = 1.05,
        maximum_crossflow_angle = 1.2)
    minimum_shape_factor > 1.0 ||
        throw(ArgumentError("minimum_shape_factor must exceed 1"))
    initial_shape_factor > minimum_shape_factor ||
        throw(ArgumentError("initial_shape_factor must exceed minimum_shape_factor"))
    separation_shape_factor > initial_shape_factor ||
        throw(ArgumentError("separation_shape_factor must exceed initial_shape_factor"))
    maximum_crossflow_angle > 0 ||
        throw(ArgumentError("maximum_crossflow_angle must be positive"))
    values = promote(initial_shape_factor, separation_shape_factor,
        minimum_shape_factor, maximum_crossflow_angle)
    return ThreeDimensionalClosure{typeof(first(values))}(values...)
end

raw"""
    ThreeDimensionalBoundaryLayerResult

Solution of [`solve_three_dimensional_boundary_layer`](@ref).

Per-panel fields carry the primary state and the derived quantities needed by
the viscous-inviscid coupling. `crossflow_angle` is the wall crossflow angle
``\beta_w``, measured from the local external streamline; it is identically
zero in a two-dimensional flow and is the quantity the earlier quasi-3D strip
march could not represent. `wall_shear` is the traction on the body, which is
no longer parallel to the edge velocity. `transpiration_velocity` is positive
along the outward body normal. `separated` flags panels where the attached
closure has been left behind, and `attached` is its complement, which is the
region the viscous loads may legitimately be integrated over.
"""
struct ThreeDimensionalBoundaryLayerResult{V, M, B, F}
    edge_speed::V
    momentum_thickness::V
    shape_factor::V
    kinetic_shape_factor::V
    crossflow_angle::V
    shear_coefficient::V
    displacement_thickness::V
    crossflow_displacement_thickness::V
    skin_friction_coefficient::V
    wall_shear::M
    transpiration_velocity::V
    separated::B
    attached::B
    force::F
    moment::F
end

# ---------------------------------------------------------------------------
# Scalar closure relations (Drela, incompressible turbulent)
# ---------------------------------------------------------------------------

raw"""
    turbulent_kinetic_shape_factor(shape_factor, momentum_reynolds)

Return the kinetic energy shape factor ``H^*=\theta^*/\theta`` from Drela's
incompressible turbulent correlation.
"""
function turbulent_kinetic_shape_factor(shape_factor, momentum_reynolds)
    reynolds = max(momentum_reynolds, 200.0)
    crossover = reynolds <= 400 ? 4.0 : 3.0 + 400.0 / reynolds
    limited = max(shape_factor, 1.0001)
    base = 1.505 + 4.0 / reynolds
    if limited < crossover
        return base +
               (0.165 - 1.6 / sqrt(reynolds)) * (crossover - limited)^1.6 / limited
    end
    logarithm = log(reynolds)
    excess = limited - crossover
    return base + excess^2 * (0.04 / limited +
            0.007 * logarithm / (excess + 4.0 / logarithm)^2)
end

raw"""
    turbulent_skin_friction(shape_factor, momentum_reynolds)

Return the streamwise skin-friction coefficient ``c_{f1}`` from Drela's
incompressible turbulent correlation.
"""
function turbulent_skin_friction(shape_factor, momentum_reynolds)
    reynolds = max(momentum_reynolds, 200.0)
    limited = max(shape_factor, 1.0001)
    return 0.3 * exp(-1.33 * limited) / log10(reynolds)^(1.74 + 0.31 * limited) +
           1.1e-4 * (tanh(4.0 - limited / 0.875) - 1.0)
end

"""
    turbulent_slip_velocity(shape_factor, kinetic_shape_factor)

Return Drela's equivalent slip velocity ``U_s``.
"""
function turbulent_slip_velocity(shape_factor, kinetic_shape_factor)
    limited = max(shape_factor, 1.0001)
    slip = 0.5 * kinetic_shape_factor * (1.0 - 4.0 * (limited - 1.0) / (3.0 * limited))
    return min(slip, 0.98)
end

raw"""
    equilibrium_shear_coefficient(shape_factor, kinetic_shape_factor, slip_velocity)

Return the equilibrium shear-stress coefficient ``C_{\tau,eq}``.
"""
function equilibrium_shear_coefficient(shape_factor, kinetic_shape_factor,
        slip_velocity)
    limited = max(shape_factor, 1.0001)
    return kinetic_shape_factor * 0.015 / (1.0 - slip_velocity) *
           (limited - 1.0)^3 / (limited^2 * limited)
end

raw"""
    layer_thickness(momentum_thickness, shape_factor)

Return the boundary-layer thickness ``\delta`` used by the shear-lag equation.
"""
function layer_thickness(momentum_thickness, shape_factor)
    limited = max(shape_factor, 1.0001)
    return momentum_thickness * (3.15 + 1.72 / (limited - 1.0) + limited)
end

# ---------------------------------------------------------------------------
# Three-dimensional closure (Mughal / Lokatt tensor structure, incompressible)
# ---------------------------------------------------------------------------

raw"""
    crossflow_parameter(shape_factor, crossflow_angle)

Return Mughal's crossflow parameter ``A_c=e^{1-H_k}\beta_w``, which sets the
magnitude of every crossflow-related thickness. It vanishes with the crossflow
angle, so the whole tensor closure collapses to the two-dimensional one.
"""
crossflow_parameter(shape_factor, crossflow_angle) =
    exp(1.0 - shape_factor) * crossflow_angle

raw"""
    thickness_tensors(momentum_thickness, shape_factor, kinetic_shape_factor,
                      crossflow)

Return the momentum-loss tensor ``\theta_{ij}``, the displacement vector
``\delta^*_i`` and the energy vector ``\theta^*_i`` in the streamwise-crossflow
basis, for incompressible flow.

Following Mughal's incompressible limit,

```math
\delta^*_1=H\theta_{11},\quad \delta^*_2=-A_c\delta^*_1,\quad
\theta_{21}=-A_c\theta_{11},\quad \theta_{12}=\theta_{21}-\delta^*_2,\quad
\theta_{22}=-A_c\theta_{12},
```

with the energy thicknesses built from ``E_{11}=H^*\theta_{11}`` by the same
recursion. Every crossflow component is proportional to ``A_c``.
"""
function thickness_tensors(momentum_thickness, shape_factor, kinetic_shape_factor,
        crossflow)
    theta_11 = momentum_thickness
    displacement_1 = shape_factor * theta_11
    displacement_2 = -crossflow * displacement_1
    theta_21 = -crossflow * theta_11
    theta_12 = theta_21 - displacement_2
    theta_22 = -crossflow * theta_12

    energy_11 = kinetic_shape_factor * theta_11
    energy_12 = theta_12 + crossflow * (2 * theta_11 - energy_11)
    energy_21 = -theta_22 - crossflow * energy_12
    energy_22 = -crossflow * (energy_21 - theta_22)

    return (momentum = @SMatrix([theta_11 theta_12; theta_21 theta_22]),
        displacement = @SVector([displacement_1, displacement_2]),
        energy = @SVector([energy_11 + energy_21, energy_12 + energy_22]))
end

"""
    rotate_to_local(tensor, cosine, sine)

Rotate a second-rank surface tensor from the streamwise-crossflow basis into
the local Cartesian basis, given the direction cosines of the edge velocity.
"""
function rotate_to_local(tensor, cosine, sine)
    rotation = @SMatrix [cosine -sine; sine cosine]
    return rotation * tensor * transpose(rotation)
end

"""
    rotate_vector_to_local(vector, cosine, sine)

Rotate a surface vector from the streamwise-crossflow basis into the local
Cartesian basis.
"""
function rotate_vector_to_local(vector, cosine, sine)
    return @SVector [cosine * vector[1] - sine * vector[2],
        sine * vector[1] + cosine * vector[2]]
end

# ---------------------------------------------------------------------------
# Surface topology for the marching scheme
# ---------------------------------------------------------------------------

raw"""
    SurfaceMarchTopology

Marching structure derived from a [`StructuredPanelGrid`](@ref).

`girth_order` lists the strips in geometric order around the section, so that
consecutive entries are physical neighbours. The importer interleaves the two
mirrored hull sides, so this ordering is recovered geometrically rather than
assumed. `tangent` and `binormal` are the local surface basis vectors
``\hat e_1`` (along the strip, downstream) and ``\hat e_2=\hat n\times\hat e_1``.
"""
struct SurfaceMarchTopology{I, T, M}
    girth_order::I
    station_count::Int
    streamwise_spacing::T
    arclength::T
    girth_spacing::T
    strip_width::T
    tangent::M
    binormal::M
end

"""
    girth_ordering(grid, reference_station)

Order the strips around the section by greedy nearest-neighbour chaining from
the strip whose reference-station centre has the largest transverse offset.
Ship sections are convex enough for this to recover the girth ordering, which
the mirrored-side interleaving in the importer otherwise hides.
"""
function girth_ordering(grid::StructuredPanelGrid, reference_station::Integer)
    mesh = grid.mesh
    strip_count = length(grid.strips)
    points = [
        begin
            strip = grid.strips[index]
            panel = strip[min(reference_station, length(strip))]
            (mesh.centers[panel, 2], mesh.centers[panel, 3])
        end
        for index in 1:strip_count
    ]
    remaining = Set(1:strip_count)
    current = argmax(point[1] for point in points)
    order = [current]
    delete!(remaining, current)
    while !isempty(remaining)
        distance(index) = hypot(points[index][1] - points[current][1],
            points[index][2] - points[current][2])
        current = argmin(distance, remaining)
        push!(order, current)
        delete!(remaining, current)
    end
    return order
end

"""
    build_march_topology(grid; reference_station = nothing)

Assemble the streamwise and girthwise metrics and the local surface basis used
by [`solve_three_dimensional_boundary_layer`](@ref). All strips must have equal
length, which the structured importers guarantee.
"""
function build_march_topology(grid::StructuredPanelGrid; reference_station = nothing)
    mesh = grid.mesh
    strip_count = length(grid.strips)
    station_count = length(first(grid.strips))
    all(length(strip) == station_count for strip in grid.strips) || throw(ArgumentError(
        "the three-dimensional march requires a structured grid with " *
        "equal-length strips"))
    station = isnothing(reference_station) ? cld(station_count, 2) : reference_station
    order = girth_ordering(grid, station)

    streamwise = zeros(station_count, strip_count)
    arclength = zeros(station_count, strip_count)
    girth = zeros(station_count, strip_count)
    width = zeros(station_count, strip_count)
    tangent = zeros(station_count, strip_count, 3)
    binormal = zeros(station_count, strip_count, 3)

    for (position, strip_index) in enumerate(order)
        strip = grid.strips[strip_index]
        for station_index in 1:station_count
            panel = strip[station_index]
            previous = strip[max(station_index - 1, 1)]
            following = strip[min(station_index + 1, station_count)]
            step = @view(mesh.centers[following, :]) .- @view(mesh.centers[previous, :])
            spacing = norm(step)
            spacing > 0 || throw(ArgumentError("coincident strip-panel centres"))
            divisor = (following == previous) ? 1 : count(
                !=(station_index), (max(station_index - 1, 1),
                    min(station_index + 1, station_count)))
            streamwise[station_index, position] = spacing / max(divisor, 1)

            normal = @view(mesh.normals[panel, :])
            direction = step ./ spacing
            projected = direction .- dot(direction, normal) .* normal
            projected_norm = norm(projected)
            first_axis = projected_norm > 1e-12 ? projected ./ projected_norm :
                         direction
            second_axis = cross(normal, first_axis)
            tangent[station_index, position, :] .= first_axis
            binormal[station_index, position, :] .= second_axis

            left = order[max(position - 1, 1)]
            right = order[min(position + 1, strip_count)]
            left_panel = grid.strips[left][station_index]
            right_panel = grid.strips[right][station_index]
            separation = norm(@view(mesh.centers[right_panel, :]) .-
                              @view(mesh.centers[left_panel, :]))
            neighbours = (left == strip_index ? 0 : 1) + (right == strip_index ? 0 : 1)
            girth[station_index, position] = separation / max(neighbours, 1)
            width[station_index, position] = mesh.areas[panel] /
                                             streamwise[station_index, position]
        end
    end

    for position in 1:strip_count
        arclength[1, position] = streamwise[1, position] / 2
        for station in 2:station_count
            arclength[station, position] = arclength[station - 1, position] +
                                           streamwise[station, position]
        end
    end

    return SurfaceMarchTopology(order, station_count, streamwise, arclength, girth,
        width, tangent, binormal)
end

# ---------------------------------------------------------------------------
# Station residuals and march
# ---------------------------------------------------------------------------

# The primary state is stored in transformed variables so that the Newton solve
# cannot leave the closure's validity range: the momentum thickness and the
# shape-factor excess are strictly positive, and the crossflow angle is bounded.
function unpack_state(unknowns, closure)
    momentum_thickness = exp(unknowns[1])
    shape_factor = closure.minimum_shape_factor + exp(unknowns[2])
    crossflow_angle = closure.maximum_crossflow_angle * tanh(unknowns[3])
    return momentum_thickness, shape_factor, crossflow_angle
end

function pack_state(momentum_thickness, shape_factor, crossflow_angle, closure)
    excess = max(shape_factor - closure.minimum_shape_factor, 1e-6)
    ratio = clamp(crossflow_angle / closure.maximum_crossflow_angle, -0.999, 0.999)
    return @SVector [log(momentum_thickness), log(excess), atanh(ratio)]
end

"""
    panel_fluxes(unknowns, edge_speed, edge_cosine, edge_sine, shear_coefficient,
                 viscosity, closure)

Return the streamwise and girthwise flux vectors of the three conservation
equations in the local Cartesian basis, together with the closure quantities
needed for the source terms.
"""
function panel_fluxes(unknowns, edge_speed, edge_cosine, edge_sine,
        shear_coefficient, viscosity, closure)
    momentum_thickness, shape_factor, crossflow_angle = unpack_state(unknowns, closure)
    momentum_reynolds = edge_speed * momentum_thickness / viscosity
    kinetic = turbulent_kinetic_shape_factor(shape_factor, momentum_reynolds)
    crossflow = crossflow_parameter(shape_factor, crossflow_angle)
    tensors = thickness_tensors(momentum_thickness, shape_factor, kinetic, crossflow)

    momentum_local = rotate_to_local(tensors.momentum, edge_cosine, edge_sine)
    displacement_local = rotate_vector_to_local(tensors.displacement, edge_cosine,
        edge_sine)
    energy_local = rotate_vector_to_local(tensors.energy, edge_cosine, edge_sine)

    friction = turbulent_skin_friction(shape_factor, momentum_reynolds)
    friction_local = rotate_vector_to_local(
        @SVector([friction, friction * tan(crossflow_angle)]), edge_cosine, edge_sine)
    slip = turbulent_slip_velocity(shape_factor, kinetic)
    dissipation = (1.0 + crossflow^2) *
                  (0.5 * friction * slip + shear_coefficient * (1.0 - slip))

    streamwise_flux = @SVector [momentum_local[1, 1], momentum_local[2, 1],
        energy_local[1]]
    girth_flux = @SVector [momentum_local[1, 2], momentum_local[2, 2], energy_local[2]]

    return (streamwise = streamwise_flux, girth = girth_flux,
        momentum = momentum_local, displacement = displacement_local,
        energy = energy_local, friction = friction_local,
        friction_magnitude = friction, dissipation = dissipation,
        kinetic = kinetic, slip = slip, crossflow = crossflow,
        momentum_thickness = momentum_thickness, shape_factor = shape_factor,
        crossflow_angle = crossflow_angle)
end

# Source terms of the three conservation equations, Mughal's incompressible
# form. The streamwise momentum source reduces to the von Karman equation when
# the crossflow vanishes, which the test suite checks directly.
function station_sources(closure_state, edge, gradients)
    speed = edge.speed
    momentum = closure_state.momentum
    displacement = closure_state.displacement
    energy = closure_state.energy

    momentum_source_x = momentum[1, 1] * 2 * gradients.speed_x / speed +
                        momentum[1, 2] * 2 * gradients.speed_y / speed +
                        displacement[1] * gradients.u_x / speed +
                        displacement[2] * gradients.u_y / speed -
                        0.5 * closure_state.friction[1]
    momentum_source_y = momentum[2, 1] * 2 * gradients.speed_x / speed +
                        momentum[2, 2] * 2 * gradients.speed_y / speed +
                        displacement[1] * gradients.v_x / speed +
                        displacement[2] * gradients.v_y / speed -
                        0.5 * closure_state.friction[2]
    energy_source = energy[1] * 3 * gradients.speed_x / speed +
                    energy[2] * 3 * gradients.speed_y / speed -
                    2 * closure_state.dissipation

    return @SVector [momentum_source_x, momentum_source_y, energy_source]
end

function station_residual(unknowns, context)
    state = panel_fluxes(unknowns, context.edge.speed, context.edge.cosine,
        context.edge.sine, context.shear_coefficient, context.viscosity,
        context.closure)
    sources = station_sources(state, context.edge, context.gradients)
    streamwise_derivative = (state.streamwise .- context.upstream_flux) ./
                            context.streamwise_spacing
    # The girthwise derivative is upwinded on the sign of the girthwise edge
    # velocity. Central differencing is ill-posed wherever the flow crosses the
    # strips, which on a ship hull is the whole bow: there the external
    # streamlines run round the stem at a large angle to the girth lines, so
    # information genuinely propagates from one strip to the next.
    girth_derivative = context.edge.sine >= 0 ?
                       (state.girth .- context.girth_flux_left) ./
                       context.girth_spacing :
                       (context.girth_flux_right .- state.girth) ./
                       context.girth_spacing
    return streamwise_derivative .+ girth_derivative .+ sources
end

raw"""
    solve_three_dimensional_boundary_layer(grid, edge_velocity, kinematic_viscosity;
                                           kwargs...)

March the three-dimensional integral boundary layer over a structured hull
surface.

Unlike the quasi-3D strip march in [`solve_quasi3d_boundary_layer`](@ref),
which advances a single scalar momentum thickness using ``|\boldsymbol u_e|``,
this solves the full momentum-loss tensor together with a crossflow momentum
equation and a kinetic energy equation,

```math
\frac{\partial\theta_{xx}}{\partial x}+\frac{\partial\theta_{xy}}{\partial y}+g_x=0,
\quad
\frac{\partial\theta_{yx}}{\partial x}+\frac{\partial\theta_{yy}}{\partial y}+g_y=0,
\quad
\frac{\partial\theta^*_x}{\partial x}+\frac{\partial\theta^*_y}{\partial y}+g_e=0,
```

in a local Cartesian surface basis. The consequences are that the wall shear is
no longer parallel to the edge velocity, that girthwise transport of momentum
between strips is represented, and that the shape factor is governed by an
energy balance rather than an entrainment correlation.

The system is space-marched from the bow. Each station is solved for all strips
simultaneously by Gauss-Seidel sweeps, so the girthwise coupling is implicit at
convergence; `station_sweeps` and `station_tolerance` control that inner
iteration. The shear-stress coefficient is advanced by Drela's lag equation
using the upstream value.

Girthwise boundary conditions are zero-gradient at both waterline ends, which
is the exact symmetry condition for the double-body flow.

!!! warning "Validity"
    This solver is validated for flow that is broadly aligned with the strips,
    where it reproduces the turbulent flat-plate laws to a few per cent. It is
    a *space-marching* scheme, and marching is only well posed while the
    external streamlines stay within the characteristic cone of the marching
    direction. On a ship forebody the streamlines wrap around the stem at up to
    about fifty degrees to the girth lines, and there the march fails. Making
    this work on a full hull needs the globally coupled upwind surface solve of
    Lokatt and Eller rather than a march. See
    `docs/INTEGRAL_BOUNDARY_LAYER.md`.

    The first three or four stations carry a startup transient, because the
    initial station is estimated from a flat-plate correlation rather than
    solved.
"""
function solve_three_dimensional_boundary_layer(grid::StructuredPanelGrid,
        edge_velocity::AbstractMatrix, kinematic_viscosity::Real;
        rho::Real = SETTINGS.rho, reference = (0, 0, 0),
        closure::ThreeDimensionalClosure = ThreeDimensionalClosure(),
        topology::Union{Nothing, SurfaceMarchTopology} = nothing,
        minimum_edge_speed = nothing, station_sweeps::Integer = 12,
        station_tolerance::Real = 1e-9, march_substeps::Integer = 8,
        stop_at_separation::Bool = true)
    mesh = grid.mesh
    size(edge_velocity) == (mesh.nfaces, 3) || throw(DimensionMismatch(
        "edge_velocity must have size (mesh.nfaces, 3)"))
    kinematic_viscosity > 0 ||
        throw(ArgumentError("kinematic_viscosity must be positive"))
    rho > 0 || throw(ArgumentError("rho must be positive"))
    station_sweeps > 0 || throw(ArgumentError("station_sweeps must be positive"))
    march_substeps > 0 || throw(ArgumentError("march_substeps must be positive"))

    march = isnothing(topology) ? build_march_topology(grid) : topology
    order = march.girth_order
    strip_count = length(order)
    station_count = march.station_count

    raw_speed = [norm(@view(edge_velocity[panel, :])) for panel in 1:mesh.nfaces]
    maximum_speed = maximum(raw_speed)
    maximum_speed > 0 ||
        throw(ArgumentError("at least one panel must have nonzero edge speed"))
    floor_speed = isnothing(minimum_edge_speed) ?
                  sqrt(eps(Float64)) * maximum_speed : minimum_edge_speed

    element_type = promote_type(eltype(edge_velocity), eltype(mesh.centers),
        typeof(kinematic_viscosity), typeof(rho), typeof(floor_speed),
        typeof(closure.initial_shape_factor))

    # Edge velocity resolved in the local surface basis, per station and strip.
    speed = zeros(element_type, station_count, strip_count)
    cosine = zeros(element_type, station_count, strip_count)
    sine = zeros(element_type, station_count, strip_count)
    panel_of = zeros(Int, station_count, strip_count)
    for (position, strip_index) in enumerate(order)
        strip = grid.strips[strip_index]
        for station in 1:station_count
            panel = strip[station]
            panel_of[station, position] = panel
            velocity = @view edge_velocity[panel, :]
            along = dot(velocity, @view(march.tangent[station, position, :]))
            across = dot(velocity, @view(march.binormal[station, position, :]))
            magnitude = max(hypot(along, across), floor_speed)
            speed[station, position] = magnitude
            cosine[station, position] = along / magnitude
            sine[station, position] = across / magnitude
        end
    end

    unknowns = zeros(element_type, 3, station_count, strip_count)
    shear = zeros(element_type, station_count, strip_count)
    separated = falses(station_count, strip_count)

    initial = _initial_station_state(march, speed, kinematic_viscosity, closure)
    for position in 1:strip_count
        unknowns[:, 1, position] .= initial[position]
        state = panel_fluxes(@view(unknowns[:, 1, position]), speed[1, position],
            cosine[1, position], sine[1, position], zero(element_type),
            kinematic_viscosity, closure)
        shear[1, position] = equilibrium_shear_coefficient(state.shape_factor,
            state.kinetic, turbulent_slip_velocity(state.shape_factor, state.kinetic))
    end

    _march_stations!(unknowns, shear, separated, grid, march, speed, cosine, sine,
        kinematic_viscosity, closure, station_sweeps, station_tolerance,
        march_substeps, stop_at_separation)

    return _assemble_result(grid, march, panel_of, unknowns, shear, separated, speed,
        cosine, sine, edge_velocity, kinematic_viscosity, closure, rho, reference)
end

raw"""
    equilibrium_shape_factor(momentum_reynolds; bracket = (1.2, 2.0))

Return the shape factor at which a zero-pressure-gradient turbulent layer is
self-similar, that is where ``2C_D=H^*c_f/2`` with the shear coefficient at its
equilibrium value.

Starting the march from an arbitrary shape factor injects a transient whose
relaxation rate scales as ``1/\theta``, which is violent at the first station
where the layer is thinnest. Starting on the equilibrium locus removes it.
"""
function equilibrium_shape_factor(momentum_reynolds; bracket = (1.2, 2.0))
    imbalance(shape) = begin
        kinetic = turbulent_kinetic_shape_factor(shape, momentum_reynolds)
        friction = turbulent_skin_friction(shape, momentum_reynolds)
        slip = turbulent_slip_velocity(shape, kinetic)
        shear = equilibrium_shear_coefficient(shape, kinetic, slip)
        dissipation = 0.5 * friction * slip + shear * (1.0 - slip)
        return 2 * dissipation - kinetic * friction / 2
    end
    low, high = bracket
    low_value = imbalance(low)
    high_value = imbalance(high)
    low_value * high_value > 0 && return 1.4
    for _ in 1:60
        middle = (low + high) / 2
        middle_value = imbalance(middle)
        if low_value * middle_value <= 0
            high = middle
        else
            low = middle
            low_value = middle_value
        end
    end
    return (low + high) / 2
end

# The march starts on the zero-pressure-gradient equilibrium locus at the
# flat-plate momentum thickness for the strip's own run length, with zero
# crossflow. A hull strip actually begins in a stagnation region rather than on
# a flat plate, so this is a starting estimate, not a boundary condition.
function _initial_station_state(march, speed, viscosity, closure)
    strip_count = size(speed, 2)
    return [
        begin
            arclength = march.arclength[1, position]
            reynolds = max(speed[1, position] * arclength / viscosity, 1.0e4)
            thickness = 0.036 * arclength * reynolds^(-0.2)
            momentum_reynolds = speed[1, position] * thickness / viscosity
            shape = equilibrium_shape_factor(momentum_reynolds)
            pack_state(thickness, shape, 0.0, closure)
        end
        for position in 1:strip_count
    ]
end

function _local_gradients(march, speed, cosine, sine, station, position, strip_count)
    streamwise_spacing = march.streamwise_spacing[station, position]
    girth_spacing = march.girth_spacing[station, position]
    left = max(position - 1, 1)
    right = min(position + 1, strip_count)

    speed_x = (speed[station, position] - speed[station - 1, position]) /
              streamwise_spacing
    upwind = sine[station, position] >= 0
    speed_y = upwind ? (speed[station, position] - speed[station, left]) / girth_spacing :
              (speed[station, right] - speed[station, position]) / girth_spacing

    u_here = speed[station, position] * cosine[station, position]
    v_here = speed[station, position] * sine[station, position]
    u_back = speed[station - 1, position] * cosine[station - 1, position]
    v_back = speed[station - 1, position] * sine[station - 1, position]
    u_left = speed[station, left] * cosine[station, left]
    u_right = speed[station, right] * cosine[station, right]
    v_left = speed[station, left] * sine[station, left]
    v_right = speed[station, right] * sine[station, right]

    return (speed_x = speed_x, speed_y = speed_y,
        u_x = (u_here - u_back) / streamwise_spacing,
        v_x = (v_here - v_back) / streamwise_spacing,
        u_y = upwind ? (u_here - u_left) / girth_spacing :
              (u_right - u_here) / girth_spacing,
        v_y = upwind ? (v_here - v_left) / girth_spacing :
              (v_right - v_here) / girth_spacing)
end

function _march_stations!(unknowns, shear, separated, grid, march, speed, cosine, sine,
        viscosity, closure, station_sweeps, station_tolerance, march_substeps,
        stop_at_separation)
    strip_count = size(unknowns, 3)
    station_count = march.station_count
    element_type = eltype(unknowns)

    for station in 2:station_count
        for position in 1:strip_count
            unknowns[:, station, position] .= @view unknowns[:, station - 1, position]
            shear[station, position] = shear[station - 1, position]
        end

        upstream_flux = [
            panel_fluxes(@view(unknowns[:, station - 1, position]),
                speed[station - 1, position], cosine[station - 1, position],
                sine[station - 1, position], shear[station - 1, position],
                viscosity, closure).streamwise
            for position in 1:strip_count
        ]

        for _ in 1:station_sweeps
            change = zero(element_type)
            girth_flux = [
                panel_fluxes(@view(unknowns[:, station, position]),
                    speed[station, position], cosine[station, position],
                    sine[station, position], shear[station, position],
                    viscosity, closure).girth
                for position in 1:strip_count
            ]

            for position in 1:strip_count
                if stop_at_separation && station > 2 &&
                   separated[station - 1, position]
                    separated[station, position] = true
                    continue
                end
                left = max(position - 1, 1)
                right = min(position + 1, strip_count)
                context = (edge = (speed = speed[station, position],
                        cosine = cosine[station, position],
                        sine = sine[station, position]),
                    gradients = _local_gradients(march, speed, cosine, sine, station,
                        position, strip_count),
                    upstream_flux = upstream_flux[position],
                    girth_flux_left = girth_flux[left],
                    girth_flux_right = girth_flux[right],
                    streamwise_spacing = march.streamwise_spacing[station, position],
                    girth_spacing = march.girth_spacing[station, position],
                    shear_coefficient = shear[station, position],
                    viscosity = viscosity, closure = closure)

                guess = SVector{3}(@view unknowns[:, station, position])
                upstream = SVector{3}(@view unknowns[:, station - 1, position])
                updated, converged = _advance_streamwise(upstream, context,
                    march_substeps)
                change = max(change, maximum(abs, updated .- guess))
                unknowns[:, station, position] .= updated

                state = panel_fluxes(updated, speed[station, position],
                    cosine[station, position], sine[station, position],
                    shear[station, position], viscosity, closure)
                separated[station, position] = !converged ||
                                               state.shape_factor >=
                                               closure.separation_shape_factor
                shear[station, position] = advance_shear_coefficient(
                    shear[station - 1, position], state,
                    march.streamwise_spacing[station, position])
            end
            change < station_tolerance && break
        end
    end
    return nothing
end

# A boundary layer leaving a stagnation region is stiff: the first station can
# more than double the momentum thickness in one panel. The station interval is
# therefore subdivided, exactly as the two-dimensional march does, with the
# girthwise fluxes and edge-velocity gradients frozen across the substeps.
function _advance_streamwise(upstream_state, context, substeps)
    state = upstream_state
    upstream_flux = context.upstream_flux
    remaining = context.streamwise_spacing
    sub_spacing = remaining / substeps
    minimum_spacing = context.streamwise_spacing / (substeps * 512)
    while remaining > 0
        step = min(sub_spacing, remaining)
        sub_context = merge(context,
            (upstream_flux = upstream_flux, streamwise_spacing = step))
        problem = NonlinearProblem(station_residual, state, sub_context)
        solution = NonlinearSolve.solve(problem, SimpleNewtonRaphson();
            abstol = 1e-12, reltol = 1e-12, maxiters = 50)
        if !all(isfinite, solution.u) || !_is_admissible(solution.u, context.closure)
            # A failed step means the substep was too coarse for the local
            # stiffness, not that the layer has separated. Halve and retry;
            # give up only once the step is far below the nominal one.
            sub_spacing /= 2
            sub_spacing < minimum_spacing && return state, false
            continue
        end
        state = solution.u
        remaining -= step
        upstream_flux = panel_fluxes(state, context.edge.speed, context.edge.cosine,
            context.edge.sine, context.shear_coefficient, context.viscosity,
            context.closure).streamwise
    end
    return state, true
end

# The transformed variables are unbounded, so a Newton step can wander into a
# region where the closure correlations are extrapolated far beyond their fit
# range. Reject those steps rather than propagate them.
function _is_admissible(unknowns, closure)
    momentum_thickness, shape_factor, _ = unpack_state(unknowns, closure)
    return momentum_thickness > 0 && isfinite(momentum_thickness) &&
           shape_factor < closure.separation_shape_factor * 4
end

raw"""
    advance_shear_coefficient(upstream_shear, state, spacing)

Advance Drela's shear-stress lag one station.

Written for ``q=\sqrt{C_\tau}`` the lag equation is a linear relaxation,

```math
\frac{dq}{ds}=\frac{1.05}{\delta}\left(q_{eq}-q\right),
```

which is integrated exactly over the step rather than by explicit Euler. That
matters: on a hull the layer thickness is two orders of magnitude smaller than
the panel length, so ``\Delta s/\delta\sim40`` and an explicit update is
violently unstable. The exact solution simply relaxes to equilibrium.

The sign convention here is the one the physics requires. The conservation form
quoted by Lokatt and Eller carries the opposite sign, which would drive the
shear stress away from equilibrium.
"""
function advance_shear_coefficient(upstream_shear, state, spacing)
    equilibrium = max(equilibrium_shear_coefficient(state.shape_factor, state.kinetic,
        state.slip), 0.0)
    thickness = max(layer_thickness(state.momentum_thickness, state.shape_factor),
        1e-12)
    current = max(upstream_shear, 1e-12)
    relaxation = exp(-1.05 * spacing / thickness)
    root = sqrt(equilibrium) + (sqrt(current) - sqrt(equilibrium)) * relaxation
    return max(root^2, 1e-12)
end

function _assemble_result(grid, march, panel_of, unknowns, shear, separated, speed,
        cosine, sine, edge_velocity, viscosity, closure, rho, reference)
    mesh = grid.mesh
    element_type = eltype(unknowns)
    strip_count = size(unknowns, 3)
    station_count = march.station_count

    edge_speed = zeros(element_type, mesh.nfaces)
    momentum_thickness = zeros(element_type, mesh.nfaces)
    shape_factor = zeros(element_type, mesh.nfaces)
    kinetic_shape_factor = zeros(element_type, mesh.nfaces)
    crossflow_angle = zeros(element_type, mesh.nfaces)
    shear_coefficient = zeros(element_type, mesh.nfaces)
    displacement = zeros(element_type, mesh.nfaces)
    crossflow_displacement = zeros(element_type, mesh.nfaces)
    friction = zeros(element_type, mesh.nfaces)
    wall_shear = zeros(element_type, mesh.nfaces, 3)
    transpiration = zeros(element_type, mesh.nfaces)
    panel_separated = falses(mesh.nfaces)

    mass_defect = zeros(element_type, station_count, strip_count)
    for position in 1:strip_count, station in 1:station_count
        panel = panel_of[station, position]
        state = panel_fluxes(@view(unknowns[:, station, position]),
            speed[station, position], cosine[station, position],
            sine[station, position], shear[station, position], viscosity, closure)

        edge_speed[panel] = speed[station, position]
        momentum_thickness[panel] = state.momentum_thickness
        shape_factor[panel] = state.shape_factor
        kinetic_shape_factor[panel] = state.kinetic
        crossflow_angle[panel] = state.crossflow_angle
        shear_coefficient[panel] = shear[station, position]
        displacement[panel] = state.displacement[1]
        crossflow_displacement[panel] = state.displacement[2]
        friction[panel] = state.friction_magnitude
        panel_separated[panel] = separated[station, position]

        # `state.friction` is already expressed in the local surface basis, so
        # it only needs projecting onto the global axes.
        scale = separated[station, position] ? zero(element_type) :
                0.5 * rho * speed[station, position]^2
        for axis in 1:3
            wall_shear[panel, axis] = scale *
                                      (state.friction[1] *
                                       march.tangent[station, position, axis] +
                                       state.friction[2] *
                                       march.binormal[station, position, axis])
        end

        mass_defect[station, position] = separated[station, position] ?
                                         zero(element_type) :
                                         march.strip_width[station, position] *
                                         speed[station, position] *
                                         state.displacement[1]
    end

    for position in 1:strip_count
        attached_count = something(findfirst(@view(separated[:, position])),
            station_count + 1) - 1
        attached_count > 0 || continue
        indices = 1:attached_count
        divergence = _mass_defect_divergence(
            @view(march.streamwise_spacing[indices, position]),
            @view(march.strip_width[indices, position]),
            @view(mass_defect[indices, position]))
        for station in indices
            transpiration[panel_of[station, position]] = divergence[station]
        end
    end

    force = zeros(element_type, 3)
    moment = zeros(element_type, 3)
    for panel in 1:mesh.nfaces
        panel_separated[panel] && continue
        panel_force = @view(wall_shear[panel, :]) .* mesh.areas[panel]
        force .+= panel_force
        moment .+= cross(@view(mesh.centers[panel, :]) .- reference, panel_force)
    end

    return ThreeDimensionalBoundaryLayerResult(edge_speed, momentum_thickness,
        shape_factor, kinetic_shape_factor, crossflow_angle, shear_coefficient,
        displacement, crossflow_displacement, friction, wall_shear, transpiration,
        panel_separated, .!panel_separated, force, moment)
end
