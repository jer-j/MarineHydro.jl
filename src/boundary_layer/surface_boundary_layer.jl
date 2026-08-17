import LinearSolve
import NonlinearSolve

using LinearAlgebra: cross, dot, lu, norm
using NonlinearSolve: NonlinearProblem, SimpleNewtonRaphson
using SparseArrays: SparseMatrixCSC, sparse

raw"""
    _newton_step(cache, jacobian, residual)

Solve ``J\,\delta=-r`` for a Newton step, reusing a `LinearSolve` cache.

The sparsity pattern is fixed by the upwind adjacency graph and does not change
as the iterate moves, so the symbolic analysis is worth doing once and reusing
for every step. Passing the algorithm through `LinearSolve` also makes it a
choice rather than a hard-wired `lu`, which matters at hull resolutions where a
sparse direct factorisation stops being the obvious answer.

Returns `nothing` if the solve fails or produces a non-finite step, which the
caller treats as a stalled iteration rather than an error.
"""
function _newton_step(cache, jacobian, residual)
    step = try
        if isnothing(cache[])
            problem = LinearSolve.LinearProblem(jacobian, residual)
            cache[] = LinearSolve.init(problem)
        else
            cache[].A = jacobian
            cache[].b = residual
        end
        copy(LinearSolve.solve!(cache[]).u)
    catch
        return nothing
    end
    all(isfinite, step) || return nothing
    return -step
end

raw"""
    SurfaceBoundaryLayerCache

Everything about a surface boundary-layer problem that does not change while the
nonlinear system is being solved: topology, metrics, the edge velocity resolved
in each panel's basis, and the least-squares velocity gradients.

The gradients are computed from the *prescribed* edge velocity and never from
the boundary-layer solution. That is deliberate. Lokatt and Eller abandoned a
fully coupled Newton method precisely because noisy velocity-gradient estimates
corrupt the off-diagonal Jacobian blocks; freezing them here keeps those blocks
out of the system entirely and leaves a Jacobian whose sparsity is exactly the
upwind cell-adjacency graph.
"""
struct SurfaceBoundaryLayerCache{TP, TM, T, C}
    topology::TP
    metrics::TM
    speed::Vector{T}
    cosine::Vector{T}
    sine::Vector{T}
    gradients::Vector{NamedTuple{
        (:speed_x, :speed_y, :u_x, :u_y, :v_x, :v_y), NTuple{6, T}}}
    edge_donor::Vector{Int}
    edge_influx::Vector{T}
    active::Vector{Bool}
    diffusion::T
    upwind_smoothing::T
    reference_speed::T
    viscosity::T
    closure::C
end

"""
    rotation_into(topology, metrics, edge, panel)

The two-by-two matrix mapping in-surface components from the panel on the other
side of `edge` into the basis of `panel`.
"""
function rotation_into(topology::SurfaceTopology, metrics::SurfaceMetrics,
        edge::Integer, panel::Integer)
    rotation = @view metrics.edge_rotation[edge, :, :]
    if topology.edge_cells[edge, 1] == panel
        return @SMatrix [rotation[1, 1] rotation[1, 2]; rotation[2, 1] rotation[2, 2]]
    end
    return @SMatrix [rotation[1, 1] rotation[2, 1]; rotation[1, 2] rotation[2, 2]]
end

"""
    edge_partner(topology, edge, panel)

The panel on the other side of `edge`, or `0` on a boundary edge.
"""
function edge_partner(topology::SurfaceTopology, edge::Integer, panel::Integer)
    topology.edge_cells[edge, 1] == panel && return topology.edge_cells[edge, 2]
    topology.edge_cells[edge, 2] == panel && return topology.edge_cells[edge, 1]
    throw(ArgumentError("panel $(panel) does not own edge $(edge)"))
end

raw"""
    build_surface_cache(mesh, edge_velocity, kinematic_viscosity; kwargs...)

Resolve the edge velocity into each panel's local basis, compute the
least-squares velocity gradients, and choose the upwind donor of every edge.

The donor of an interior edge is the panel the flow comes *from*, decided on the
sign of the inverse-distance-weighted edge velocity resolved on the edge normal.
This is the first-order upwinding that makes the scheme stable irrespective of
how the mesh lines happen to lie relative to the flow — which is exactly what a
strip march cannot provide.
"""
function build_surface_cache(mesh::Mesh, edge_velocity::AbstractMatrix,
        kinematic_viscosity::Real; closure = ThreeDimensionalClosure(),
        topology = nothing, metrics = nothing, minimum_edge_speed = nothing,
        weld_tolerance = nothing, gradient_bound::Real = 1.0, active = nothing,
        diffusion::Real = 0.0, upwind_smoothing::Real = 0.0)
    size(edge_velocity) == (mesh.nfaces, 3) || throw(DimensionMismatch(
        "edge_velocity must have size (mesh.nfaces, 3)"))
    kinematic_viscosity > 0 ||
        throw(ArgumentError("kinematic_viscosity must be positive"))
    gradient_bound > 0 || throw(ArgumentError("gradient_bound must be positive"))

    surface = isnothing(topology) ?
              build_surface_topology(mesh; weld_tolerance) : topology
    geometry = isnothing(metrics) ? build_surface_metrics(mesh, surface) : metrics

    element_type = promote_type(eltype(edge_velocity), eltype(mesh.vertices),
        typeof(kinematic_viscosity))
    raw_speed = [norm(@view(edge_velocity[panel, :])) for panel in 1:mesh.nfaces]
    largest = maximum(raw_speed)
    largest > 0 ||
        throw(ArgumentError("at least one panel must have nonzero edge speed"))
    floor_speed = isnothing(minimum_edge_speed) ?
                  sqrt(eps(Float64)) * largest : minimum_edge_speed

    speed = Vector{element_type}(undef, mesh.nfaces)
    cosine = Vector{element_type}(undef, mesh.nfaces)
    sine = Vector{element_type}(undef, mesh.nfaces)
    for panel in 1:mesh.nfaces
        velocity = @view edge_velocity[panel, :]
        along = dot(velocity, row3(geometry.tangent, panel))
        across = dot(velocity, row3(geometry.binormal, panel))
        magnitude = max(hypot(along, across), floor_speed)
        speed[panel] = magnitude
        cosine[panel] = along / magnitude
        sine[panel] = across / magnitude
    end

    solved = if isnothing(active)
        fill(true, mesh.nfaces)
    else
        length(active) == mesh.nfaces || throw(DimensionMismatch(
            "active must have one entry per panel"))
        any(active) || throw(ArgumentError("at least one panel must be active"))
        collect(Bool, active)
    end

    gradients = _least_squares_gradients(mesh, surface, geometry, speed, cosine, sine,
        gradient_bound, solved)
    donor, influx = _upwind_donors(mesh, surface, geometry, speed, cosine, sine, solved)

    return SurfaceBoundaryLayerCache(surface, geometry, speed, cosine, sine,
        gradients, donor, influx, solved, convert(element_type, diffusion),
        convert(element_type, upwind_smoothing), convert(element_type, largest), convert(element_type, kinematic_viscosity),
        closure)
end

# Weighted linear least squares for the in-surface gradient of the edge velocity
# and its magnitude. The moment matrix depends only on geometry, so the whole
# stencil could be cached; it is rebuilt here because the metrics carry
# derivatives and the cost is negligible next to the Newton solve.
function _least_squares_gradients(mesh, topology, metrics, speed::AbstractVector{T},
        cosine, sine, bound, active) where {T}
    # `T` is a type PARAMETER, not a value argument. Passing the element type as
    # a value leaves every `SMatrix{2,2,T}` in the loop below to be built at run
    # time, which heap-allocates each one; on KVLCC2 that alone accounted for
    # most of a hundred thousand allocations per cache build.
    element_type = T
    gradient_type = NamedTuple{
        (:speed_x, :speed_y, :u_x, :u_y, :v_x, :v_y), NTuple{6, T}}
    gradients = Vector{gradient_type}(undef, mesh.nfaces)
    for panel in 1:mesh.nfaces
        moment = zero(SMatrix{2, 2, element_type})
        speed_right = zero(SVector{2, element_type})
        u_right = zero(SVector{2, element_type})
        v_right = zero(SVector{2, element_type})
        own_u = speed[panel] * cosine[panel]
        own_v = speed[panel] * sine[panel]
        for side in 1:4
            edge = topology.cell_edges[panel, side]
            edge == 0 && continue
            topology.edge_kind[edge] === :interior || continue
            other = edge_partner(topology, edge, panel)
            # A neighbour outside the solved region is outside the domain, so it
            # must not enter the stencil any more than the surface's own edge
            # would.
            active[other] || continue
            offset3 = row3(mesh.centers, other) - row3(mesh.centers, panel)
            offset = @SVector [dot(offset3, row3(metrics.tangent, panel)),
                dot(offset3, row3(metrics.binormal, panel))]
            distance = hypot(offset[1], offset[2])
            distance > 0 || continue
            weight = one(element_type) / distance^2
            rotation = rotation_into(topology, metrics, edge, panel)
            neighbour = rotation *
                        @SVector([speed[other] * cosine[other],
                speed[other] * sine[other]])
            moment = moment + weight * (offset * transpose(offset))
            speed_right = speed_right + weight * (speed[other] - speed[panel]) * offset
            u_right = u_right + weight * (neighbour[1] - own_u) * offset
            v_right = v_right + weight * (neighbour[2] - own_v) * offset
        end
        determinant = moment[1, 1] * moment[2, 2] - moment[1, 2] * moment[2, 1]
        if abs(determinant) < eps(element_type) * max(one(element_type),
            abs(moment[1, 1] * moment[2, 2]))
            gradients[panel] = gradient_type((zero(element_type), zero(element_type),
                zero(element_type), zero(element_type), zero(element_type),
                zero(element_type)))
            continue
        end
        inverse = SMatrix{2, 2, element_type}(moment[2, 2] / determinant,
            -moment[2, 1] / determinant, -moment[1, 2] / determinant,
            moment[1, 1] / determinant)
        speed_gradient = inverse * speed_right
        u_gradient = inverse * u_right
        v_gradient = inverse * v_right
        rescale = _gradient_rescaling(speed_gradient, u_gradient, v_gradient,
            sqrt(mesh.areas[panel]), speed[panel], bound)
        gradients[panel] = gradient_type((rescale * speed_gradient[1],
            rescale * speed_gradient[2], rescale * u_gradient[1],
            rescale * u_gradient[2], rescale * v_gradient[1],
            rescale * v_gradient[2]))
    end
    return gradients
end

raw"""
    _gradient_rescaling(speed_gradient, u_gradient, v_gradient, size, speed, bound)

Factor shrinking a cell's velocity gradients until the largest of them changes
the edge velocity by no more than `bound` times itself across one cell.

This is Lokatt's bounding of the velocity-gradient terms, and without it the
solve does not run on a hull at all. Around the stem the local basis swings
sharply from panel to panel across cells that are slivers, so differencing the
rotated velocity components over short centroid separations produces gradients
two orders of magnitude larger than anywhere else — on the KVLCC2 forebody
``\partial u/\partial x`` reaches 129 against 0.04 on an ordinary panel. Those
feed straight into `station_sources`, whose streamwise term then demands more
momentum than the cell's inflow can supply, and the cell's momentum equation
loses its positive-thickness root entirely. Newton is then asked to solve a
system a handful of whose rows have no admissible solution, and it stalls with
the residual concentrated on a thousandth of the wetted area.

Rescaling uniformly rather than clipping each component keeps the direction of
the gradient, which is what carries the crossflow information. Lokatt notes the
same device perturbs the initial shape-factor development near the attachment
line; that is the price of getting a solve at all there.
"""
function _gradient_rescaling(speed_gradient, u_gradient, v_gradient, size, speed,
        bound)
    # Frobenius norms, not the largest component. The local basis is fixed by
    # the panel's own vertex order and has no physical orientation, so a measure
    # that depends on which way it happens to point makes the rescaling — and
    # with it the answer — a function of the mesh's winding. On a hull whose two
    # halves are mirrored that shows up immediately as an asymmetric solution to
    # a symmetric problem.
    velocity = sqrt(u_gradient[1]^2 + u_gradient[2]^2 + v_gradient[1]^2 +
                    v_gradient[2]^2)
    magnitude = hypot(speed_gradient[1], speed_gradient[2])
    scaled = max(velocity, magnitude) * size / max(speed, eps(typeof(speed)))
    return scaled > bound ? bound / scaled : one(scaled)
end

function _upwind_donors(mesh, topology, metrics, speed, cosine, sine, active)
    donor = zeros(Int, topology.nedges)
    influx = zeros(eltype(speed), topology.nedges)
    for edge in 1:topology.nedges
        left = topology.edge_cells[edge, 1]
        right = topology.edge_cells[edge, 2]
        # A cell outside the solved region acts exactly like the outside of the
        # surface: the edge between it and a solved cell becomes a boundary of
        # the domain, through which the defect can leave. That is what gives the
        # convergence region at the stern somewhere to discharge into.
        left_solved = left != 0 && active[left]
        right_solved = right != 0 && active[right]
        (left_solved || right_solved) || continue
        owner = right_solved && !left_solved ? right : left
        if !(left_solved && right_solved)
            # On a boundary edge, outflow extrapolates from the owning cell and
            # inflow is marked with a zero donor so that a prescribed state is
            # supplied instead. Without that the flux divergence of a leading
            # cell would cancel by the geometric conservation law and the
            # momentum equation would have no solution.
            outward = _outward_speed(metrics, edge, owner, topology, speed, cosine,
                sine)
            donor[edge] = outward >= 0 ? owner : 0
            influx[edge] = abs(outward)
            continue
        end
        normal = (metrics.edge_area_vector[edge, 1, 1],
            metrics.edge_area_vector[edge, 1, 2])
        left_velocity = (speed[left] * cosine[left], speed[left] * sine[left])
        rotation = rotation_into(topology, metrics, edge, left)
        right_velocity = rotation *
                         @SVector([speed[right] * cosine[right],
            speed[right] * sine[right]])
        left_weight = one(eltype(speed)) /
                      max(norm(row3(metrics.edge_midpoint, edge) -
                               row3(mesh.centers, left)), eps(Float64))
        right_weight = one(eltype(speed)) /
                       max(norm(row3(metrics.edge_midpoint, edge) -
                                row3(mesh.centers, right)), eps(Float64))
        blended = ((left_weight * left_velocity[1] + right_weight * right_velocity[1]),
            (left_weight * left_velocity[2] + right_weight * right_velocity[2])) ./
                  (left_weight + right_weight)
        outward = blended[1] * normal[1] + blended[2] * normal[2]
        donor[edge] = outward >= 0 ? left : right
        influx[edge] = abs(outward)
    end
    return donor, influx
end

# Edge velocity of `owner` resolved on the outward edge normal it sees.
function _outward_speed(metrics, edge, owner, topology, speed, cosine, sine)
    area_vector = edge_area_vector_for(topology, metrics, edge, owner)
    return speed[owner] * (cosine[owner] * area_vector[1] +
            sine[owner] * area_vector[2])
end

raw"""
    attached_flow_domain(mesh, edge_velocity, kinematic_viscosity; kwargs...)

The largest part of the surface on which the boundary-layer equations have an
outlet, as a per-panel mask for the `active` keyword.

A closed double body has none globally: every finite acyclic upwind graph has a
sink, and on KVLCC2 that sink is a cell at the stern which takes flux from all
four sides and donates to none. Its own state enters its own residual only
through the source term, of order ``c_f/2``, while it must absorb the whole
defect arriving from upstream — its equation has no solution, and a global line
search lets it stall the entire hull.

That cell is not a discretisation artifact. The surface streamlines genuinely
converge there and the fluid leaves the surface, which is three-dimensional
separation. Tanaka is explicit that first-order integral methods do not apply
near the stern end of a full-form ship, where the boundary layer becomes thick
and bilge vortices form by exactly this mechanism, and that the loads there need
a separately modelled longitudinal vortex superposed on the layer. So the right
response is to stop the integral calculation where its assumptions fail, not to
force a solution out of it.

Cells without an outlet are removed and the donors recomputed, repeatedly. That
terminates quickly, because removing a sink turns its edges into domain
boundaries and so hands the cells that fed it an outflow of their own. What is
left is a domain on which every cell can discharge.

Returns the mask and the fraction of wetted area retained. A retained fraction
far below one says the flow separates over much of the hull and the integral
result should not be trusted as a whole-body answer.
"""
function attached_flow_domain(mesh::Mesh, edge_velocity::AbstractMatrix,
        kinematic_viscosity::Real; topology = nothing, metrics = nothing,
        minimum_retained::Real = 0.5, closure = ThreeDimensionalClosure(),
        curvature_limit::Real = 0.5, kwargs...)
    surface = isnothing(topology) ? build_surface_topology(mesh) : topology
    geometry = isnothing(metrics) ? build_surface_metrics(mesh, surface) : metrics
    active = fill(true, mesh.nfaces)

    # Thin-layer validity first. Everything the integral equations assume rests
    # on the layer being thin against the local radius of curvature, and where
    # that fails the equations are the wrong ones rather than merely hard to
    # solve. Tanaka calls the ordinary formulation "first-level boundary layer
    # theory", notes that local curvatures are not in it, and that where the
    # radius is small "the basic governing equations become different from the
    # ordinary ones" — a thick boundary layer problem compared with the local
    # radius of curvature.
    if curvature_limit > 0
        thickness = _estimated_thickness(mesh, surface, geometry, edge_velocity,
            kinematic_viscosity, closure)
        for panel in 1:mesh.nfaces
            thickness[panel] * surface_curvature(mesh, surface, panel) >
            curvature_limit || continue
            active[panel] = false
        end
    end
    for _ in 1:(mesh.nfaces)
        cache = build_surface_cache(mesh, edge_velocity, kinematic_viscosity;
            topology = surface, metrics = geometry, active, closure, kwargs...)
        removed = 0
        for panel in 1:mesh.nfaces
            active[panel] || continue
            outlet = false
            for side in 1:4
                edge = surface.cell_edges[panel, side]
                edge == 0 && continue
                if cache.edge_donor[edge] == panel
                    outlet = true
                    break
                end
            end
            outlet && continue
            active[panel] = false
            removed += 1
        end
        removed == 0 && break
    end

    retained = sum(mesh.areas[active]) / sum(mesh.areas)
    retained >= minimum_retained || @warn "attached_flow_domain retained only " *
        "$(round(100 * retained, digits = 1))% of the wetted area; the integral " *
        "boundary layer separates over most of this hull"
    return active, retained
end

raw"""
    surface_curvature(mesh, topology, panel)

Largest rate of change of the surface normal across a panel's interior edges,
an estimate of the local principal curvature ``\kappa=1/R``.

Taken as a maximum rather than a mean because the thin-layer assumption fails
as soon as it fails in *any* direction: a bilge is sharply curved girthwise
while being nearly straight along the hull.
"""
function surface_curvature(mesh::Mesh, topology::SurfaceTopology, panel::Integer)
    largest = zero(eltype(mesh.centers))
    for side in 1:4
        edge = topology.cell_edges[panel, side]
        edge == 0 && continue
        topology.edge_kind[edge] === :interior || continue
        other = edge_partner(topology, edge, panel)
        other == 0 && continue
        separation = norm(@view(mesh.centers[other, :]) .-
                          @view(mesh.centers[panel, :]))
        separation > 0 || continue
        turn = norm(row3(mesh.normals, other) - row3(mesh.normals, panel))
        largest = max(largest, turn / separation)
    end
    return largest
end

# Boundary-layer thickness estimated from the flat-plate correlation at each
# panel's arclength, so the validity test needs no solve and cannot depend on
# one that failed.
function _estimated_thickness(mesh, topology, metrics, edge_velocity,
        kinematic_viscosity, closure)
    cache = build_surface_cache(mesh, edge_velocity, kinematic_viscosity; closure,
        topology, metrics)
    order = flow_ordering(cache)
    states = initial_states(mesh, cache, order, inflow_states(mesh, cache))
    thickness = zeros(eltype(cache.speed), mesh.nfaces)
    for panel in 1:mesh.nfaces
        momentum, shape, _ = unpack_state(SVector{3}(@view states[:, panel]), closure)
        thickness[panel] = layer_thickness(momentum, shape)
    end
    return thickness
end

raw"""
    edge_flux_contribution(donor_state, edge, panel, cache)

Flux of the three conserved quantities through `edge` into `panel`, already
weighted by the edge area vector.

The momentum defect is a second-rank surface tensor with *both* indices in the
surface, so transporting it from the donor's basis into the receiving panel's
requires ``R\Theta R^{\mathsf T}``, not ``R\Theta``. The energy defect is a
vector and takes ``RE``. Getting the tensor rank wrong is silent and produces a
plausible but wrong answer, which is why the test suite checks invariance under
a rigid rotation of the whole problem.
"""
function edge_flux_contribution(donor_state, edge::Integer, panel::Integer,
        cache::SurfaceBoundaryLayerCache, shear, opposite_state = nothing)
    donor = cache.edge_donor[edge]
    # An inflow boundary edge carries the prescribed state, evaluated in the
    # receiving panel's own frame.
    source_panel = donor == 0 ? panel : donor
    momentum, energy = _edge_defects(donor_state, source_panel, edge, panel, cache,
        shear)
    weight = _upwind_weight(cache, edge, panel)
    if weight < 1 && !isnothing(opposite_state)
        # Blend in the other side. See `_upwind_weight`.
        other = edge_partner(cache.topology, edge, source_panel)
        other = other == 0 ? panel : other
        opposite_momentum, opposite_energy = _edge_defects(opposite_state, other, edge,
            panel, cache, shear)
        momentum = weight * momentum + (1 - weight) * opposite_momentum
        energy = weight * energy + (1 - weight) * opposite_energy
    end
    area_vector = edge_area_vector_for(cache.topology, cache.metrics, edge, panel)
    normal = @SVector [area_vector[1], area_vector[2]]
    transported = momentum * normal
    return @SVector [transported[1], transported[2], dot(energy, normal)]
end

# Momentum and energy defect of `source_panel`, rotated into `panel`'s basis.
function _edge_defects(state_vector, source_panel, edge, panel, cache, shear)
    state = panel_fluxes(state_vector, cache.speed[source_panel],
        cache.cosine[source_panel], cache.sine[source_panel], shear[source_panel],
        cache.viscosity, cache.closure)
    momentum = state.momentum
    energy = state.energy
    if source_panel != panel
        rotation = rotation_into(cache.topology, cache.metrics, edge, panel)
        momentum = rotation * momentum * transpose(rotation)
        energy = rotation * energy
    end
    return momentum, energy
end

raw"""
    _upwind_weight(cache, edge, panel)

Share of an edge flux taken from the upwind side, going smoothly to one half as
the flow becomes parallel to the edge.

```math
w=\tfrac12\left[1+\tanh\!\left(\frac{|q\cdot m|}{\epsilon_u\,q\,\ell}\right)\right]
```

A hard upwind switch is a discrete function of the external velocity, so a
Newton method cannot see that moving the external flow would flip it. That is
fatal in a strongly coupled solve at precisely the place it matters: a
convergence line is where ``q\cdot m\to0`` and the switch is degenerate. With
`upwind_smoothing` at zero this returns one and the scheme is the plain
first-order upwind it was; positive values recover a central average exactly
where the flow runs along an edge, which is the same weakly-elliptic behaviour
Drela obtains from added diffusion, and which is what lets information cross a
convergence line.

It also removes a discrete branch from the map, which is what mesh sensitivities
need.
"""
function _upwind_weight(cache::SurfaceBoundaryLayerCache, edge::Integer,
        panel::Integer)
    smoothing = cache.upwind_smoothing
    smoothing > 0 || return one(smoothing)
    area_vector = edge_area_vector_for(cache.topology, cache.metrics, edge, panel)
    length_scale = max(hypot(area_vector[1], area_vector[2]), eps(Float64))
    flux = abs(cache.speed[panel] * (cache.cosine[panel] * area_vector[1] +
                cache.sine[panel] * area_vector[2]))
    return (one(smoothing) +
            tanh(flux / (smoothing * cache.speed[panel] * length_scale))) / 2
end

raw"""
    cell_residual(own_state, donor_states, panel, mesh, cache)

Residual of the three conservation equations on `panel`,

```math
R_i=\frac{1}{A_i}\sum_{e\in\partial i}\ell_e\,\Phi_e^{(i)}
+g\left(U_i,\nabla q_i\right),
```

with `donor_states` giving the primary state of each side's upwind donor in the
order of `cell_edges`.
"""
function cell_residual(own_state, donor_states, panel::Integer, mesh::Mesh,
        cache::SurfaceBoundaryLayerCache, shear, neighbour_states = nothing)
    flux_sum = zero(SVector{3, eltype(own_state)})
    for side in 1:4
        edge = cache.topology.cell_edges[panel, side]
        edge == 0 && continue
        opposite = isnothing(neighbour_states) ? nothing : neighbour_states[side]
        flux_sum = flux_sum +
                   edge_flux_contribution(donor_states[side], edge, panel, cache,
            shear, opposite)
    end
    if !isnothing(neighbour_states)
        flux_sum = flux_sum + _diffusive_flux(own_state, neighbour_states, panel,
            mesh, cache, shear)
    end
    state = panel_fluxes(own_state, cache.speed[panel], cache.cosine[panel],
        cache.sine[panel], shear[panel], cache.viscosity, cache.closure)
    sources = station_sources(state, (speed = cache.speed[panel],), cache.gradients[panel])
    return flux_sum ./ mesh.areas[panel] .+ sources
end

raw"""
    _diffusive_flux(own_state, neighbour_states, panel, mesh, cache, shear)

Drela's added conservative numerical diffusion, as an edge flux.

```math
\nabla\cdot\left[V_\epsilon\,\bar{\bar h}\,\nabla M\right],
\qquad V_\epsilon=\epsilon\max_j q_j
```

Applied to the same defect fluxes the convective term carries, contracted on
the same edge area vector, so it is antisymmetric between the two cells sharing
an edge and therefore conservative: it may smear a gradient but cannot change
the total momentum or energy defect leaving the domain.

This is the term that lets a *convergence line* be represented. Where the
external streamlines converge, the characteristics of the hyperbolic system
converge too, and Drela notes such lines are captured "with no special treatment
being required", analogous to a captured shock. Pure upwinding cannot do this:
at a cell where every edge is an inflow there is no downstream direction to
smear into, its equation has no solution, and the solve stalls on it. The
diffusion term renders the system weakly elliptic, which gives information a way
across the line.

The cost is the block-lower-triangular structure, which this term destroys — it
couples each cell to its neighbours in both directions. That is the honest trade
between a purely hyperbolic scheme that cannot represent separation and a weakly
elliptic one that can.
"""
function _diffusive_flux(own_state, neighbour_states, panel::Integer, mesh::Mesh,
        cache::SurfaceBoundaryLayerCache, shear)
    total = zero(SVector{3, eltype(own_state)})
    cache.diffusion > 0 || return total
    own = panel_fluxes(own_state, cache.speed[panel], cache.cosine[panel],
        cache.sine[panel], shear[panel], cache.viscosity, cache.closure)
    strength = cache.diffusion * cache.reference_speed
    for side in 1:4
        edge = cache.topology.cell_edges[panel, side]
        edge == 0 && continue
        cache.topology.edge_kind[edge] === :interior || continue
        other = edge_partner(topology_of(cache), edge, panel)
        (other == 0 || !cache.active[other]) && continue
        neighbour = panel_fluxes(neighbour_states[side], cache.speed[other],
            cache.cosine[other], cache.sine[other], shear[other], cache.viscosity,
            cache.closure)
        rotation = rotation_into(cache.topology, cache.metrics, edge, panel)
        momentum = rotation * neighbour.momentum * transpose(rotation) - own.momentum
        energy = rotation * neighbour.energy - own.energy
        area_vector = edge_area_vector_for(cache.topology, cache.metrics, edge, panel)
        normal = @SVector [area_vector[1], area_vector[2]]
        separation = max(norm(row3(mesh.centers, other) - row3(mesh.centers, panel)),
            eps(Float64))
        # The grid length matrix reduces here to the cell size along the edge
        # normal, which for a quadrilateral is the centroid separation.
        weight = strength * hypot(area_vector[1], area_vector[2]) / separation
        transported = momentum * normal
        total = total + weight *
                        @SVector [transported[1], transported[2], dot(energy, normal)]
    end
    return total
end

# Both the smoothed upwind blend and the diffusion term read the panel across
# each side, not just the donor.
_needs_neighbours(cache::SurfaceBoundaryLayerCache) =
    cache.diffusion > 0 || cache.upwind_smoothing > 0

topology_of(cache::SurfaceBoundaryLayerCache) = cache.topology

# Primary states of the donor of each of a panel's sides, taken from the global
# state array. A boundary edge donates the panel's own state, which makes an
# outflow edge a first-order extrapolation.
function _donor_states(states, panel, cache, inflow)
    return ntuple(4) do side
        edge = cache.topology.cell_edges[panel, side]
        edge == 0 && return zero(SVector{3, eltype(states)})
        donor = cache.edge_donor[edge]
        donor == 0 && return inflow[edge]
        return SVector{3}(@view states[:, donor])
    end
end

# States of the panel across each side, irrespective of upwind direction. The
# diffusion term is symmetric, so it needs both neighbours, not just the donor.
function _neighbour_states(states, panel, cache)
    return ntuple(4) do side
        edge = cache.topology.cell_edges[panel, side]
        edge == 0 && return zero(SVector{3, eltype(states)})
        cache.topology.edge_kind[edge] === :interior ||
            return SVector{3}(@view states[:, panel])
        other = edge_partner(cache.topology, edge, panel)
        (other == 0 || !cache.active[other]) &&
            return SVector{3}(@view states[:, panel])
        return SVector{3}(@view states[:, other])
    end
end

raw"""
    _inflow_strength(cache, panel, edge)

Volume of flow entering `panel` through `edge`, and zero when the edge is an
outflow or the flow merely runs along it.

Read straight off the donor chosen by `_upwind_donors`, so that the two can
never disagree. They must not: the donor is decided on the inverse-distance
*blended* edge velocity, and around the stem the two neighbouring cells point
differently enough that a panel's own velocity gives the opposite sign. When
this was computed independently from the panel's own velocity, `flow_ordering`
was therefore not a topological order of the graph the fluxes actually use, and
`initial_states` propagated arclength backwards into the stem — giving cells at
the very start of the hull a momentum thickness four times the mid-body value.

A side the flow runs parallel to carries a donor by the `>= 0` tie-break but no
flux, and its influx is zero to rounding, which is what keeps a mesh aligned
with the stream from chaining its rows together.
"""
function _inflow_strength(cache::SurfaceBoundaryLayerCache, panel::Integer,
        edge::Integer)
    cache.edge_donor[edge] == panel && return zero(eltype(cache.edge_influx))
    influx = cache.edge_influx[edge]
    area_vector = edge_area_vector_for(cache.topology, cache.metrics, edge, panel)
    threshold = sqrt(eps(Float64)) * cache.speed[panel] *
                hypot(area_vector[1], area_vector[2])
    return influx > threshold ? influx : zero(influx)
end

raw"""
    inflow_states(mesh, cache; run_length = nothing)

Prescribed primary state on every inflow boundary edge.

A boundary edge that the flow enters through needs a genuine condition: with
zero-gradient there the flux divergence of a leading cell cancels by the
geometric conservation law and the streamwise momentum equation has no
solution. The default is the turbulent flat-plate estimate at half a cell,
sitting on the zero-pressure-gradient equilibrium locus with no crossflow —
the same starting state the earlier strip march used at its first station.
"""
function inflow_states(mesh::Mesh, cache::SurfaceBoundaryLayerCache;
        run_length = nothing)
    element_type = eltype(cache.speed)
    states = Vector{SVector{3, element_type}}(undef, cache.topology.nedges)
    default = zero(SVector{3, element_type})
    fill!(states, default)
    for edge in 1:cache.topology.nedges
        cache.edge_donor[edge] == 0 || continue
        panel = cache.topology.edge_cells[edge, 1]
        distance = isnothing(run_length) ? sqrt(mesh.areas[panel]) / 2 : run_length
        reynolds = max(cache.speed[panel] * distance / cache.viscosity, 1.0e4)
        thickness = 0.036 * distance * reynolds^(-0.2)
        momentum_reynolds = cache.speed[panel] * thickness / cache.viscosity
        shape = equilibrium_shape_factor(momentum_reynolds)
        states[edge] = pack_state(thickness, shape, zero(element_type), cache.closure)
    end
    return states
end

raw"""
    initial_states(mesh, cache, order, inflow)

Starting guess for the global solve.

An arclength is propagated along the flow ordering and the flat-plate
correlation at that arclength gives the momentum thickness, with the shape
factor on the equilibrium locus and no crossflow. This is only a starting
estimate; it is not a boundary condition, and the converged answer does not
depend on it.

Each panel takes the *inflow-weighted mean* of its upwind neighbours' arclength
plus the centroid separation resolved along the flow. Weighting by
``\max(-q\cdot m,0)`` is what keeps a side the flow merely runs parallel to out
of the sum: on a plate meshed square to a uniform stream those sides carry
``q\cdot m=0`` exactly, and a plain maximum over neighbours would add the
transverse cell size once per row, doubling the arclength halfway across.
"""
function initial_states(mesh::Mesh, cache::SurfaceBoundaryLayerCache,
        order::AbstractVector{<:Integer}, inflow)
    element_type = eltype(cache.speed)
    arclength = zeros(element_type, mesh.nfaces)
    for panel in order
        weighted = zero(element_type)
        total = zero(element_type)
        for side in 1:4
            edge = cache.topology.cell_edges[panel, side]
            edge == 0 && continue
            donor = cache.edge_donor[edge]
            (donor == 0 || donor == panel) && continue
            weight = _inflow_strength(cache, panel, edge)
            weight > 0 || continue
            offset = row3(mesh.centers, panel) - row3(mesh.centers, donor)
            step = abs(dot(offset, row3(cache.metrics.tangent, panel)) *
                       cache.cosine[panel] +
                       dot(offset, row3(cache.metrics.binormal, panel)) *
                       cache.sine[panel])
            weighted += weight * (arclength[donor] + step)
            total += weight
        end
        arclength[panel] = total > 0 ? weighted / total : sqrt(mesh.areas[panel]) / 2
    end

    states = zeros(element_type, 3, mesh.nfaces)
    for panel in 1:mesh.nfaces
        reynolds = max(cache.speed[panel] * arclength[panel] / cache.viscosity, 1.0e4)
        thickness = 0.036 * arclength[panel] * reynolds^(-0.2)
        momentum_reynolds = cache.speed[panel] * thickness / cache.viscosity
        shape = equilibrium_shape_factor(momentum_reynolds)
        states[:, panel] .= pack_state(thickness, shape, zero(element_type),
            cache.closure)
    end
    return states
end

raw"""
    flow_ordering(cache)

Panels sorted so that every panel follows its upwind neighbours wherever the
upwind graph allows it, by Kahn's algorithm on that graph.

This is the piece that replaces the strip march. Sweeping in this order makes a
nonlinear Gauss-Seidel pass behave like a march, but along the *actual* flow
rather than along mesh lines, so it stays well posed where streamlines cross the
girth lines at sixty degrees. Cycles — which recirculation produces — are broken
by emitting the remaining panels in their existing order; the global Newton
solve that follows does not depend on the ordering being exact.
"""
function flow_ordering(cache::SurfaceBoundaryLayerCache)
    topology = cache.topology
    nfaces = size(topology.cell_edges, 1)
    upstream_count = zeros(Int, nfaces)
    downstream = [Int[] for _ in 1:nfaces]
    for edge in 1:topology.nedges
        topology.edge_kind[edge] === :interior || continue
        donor = cache.edge_donor[edge]
        # An interior edge onto a deactivated panel is a domain boundary and
        # carries a prescribed inflow, so it orders nothing.
        donor == 0 && continue
        receiver = edge_partner(topology, edge, donor)
        (receiver == 0 || !cache.active[receiver]) && continue
        # Every edge whose donor is the other cell, whether or not it carries
        # flux. This must be the residual's dependency graph exactly: a
        # zero-flux edge still supplies a donor state, and its flux still has a
        # nonzero derivative, so leaving it out gave an order in which the
        # permuted Jacobian was not triangular and a sweep that was not a march.
        push!(downstream[donor], receiver)
        upstream_count[receiver] += 1
    end
    ready = [panel for panel in 1:nfaces
             if upstream_count[panel] == 0 && cache.active[panel]]
    order = Int[]
    sizehint!(order, nfaces)
    while !isempty(ready)
        panel = pop!(ready)
        push!(order, panel)
        for receiver in downstream[panel]
            upstream_count[receiver] -= 1
            upstream_count[receiver] == 0 && push!(ready, receiver)
        end
    end
    if length(order) < nfaces
        visited = falses(nfaces)
        visited[order] .= true
        append!(order,
            (panel for panel in 1:nfaces if !visited[panel] && cache.active[panel]))
    end
    return order
end

raw"""
    SurfaceBoundaryLayerSolution

Converged primary states and the diagnostics of the solve.
"""
struct SurfaceBoundaryLayerSolution{T}
    states::Matrix{T}
    shear_coefficient::Vector{T}
    residual_norm::T
    sweeps::Int
    newton_steps::Int
    converged::Bool
end

# Local three-by-three solve for one panel with its donors frozen. Used by the
# Gauss-Seidel globalisation, not by the Newton phase.
raw"""
    _march_panel!(states, panel, mesh, cache, inflow, shear, step_limit)

One damped Newton step on a single cell, judged on that cell's own residual.

Because the residual is first-order upwind and the flow ordering is a
topological order of its dependency graph, the Jacobian permuted into that order
is *exactly* block lower triangular. Sweeping it cell by cell is therefore a
march along the true flow rather than along mesh lines, and — the point here —
each cell's acceptance test is local and independent of every other.

That independence is what makes a closed body solvable at all. A double body
has no outlet: every finite directed acyclic graph has a sink, so some cell
receives flux from all sides and donates to none. Its own state then enters its
own residual only through the source term, of order ``c_f/2``, while it must
absorb the whole defect arriving from upstream. Its equation has no solution.
Under a *global* line search that one cell dominates the residual norm and
blocks the step for every other cell on the hull; under this one it simply fails
to improve, alone, and the rest of the surface converges around it.

Taking one damped step rather than converging the cell outright matters too: the
local three-by-three problems are stiff enough that a plain Newton solve on them
hits its iteration cap on half the cells of a flat plate.
"""
function _march_panel!(states, panel, mesh, cache, inflow, shear, step_limit)
    donors = _donor_states(states, panel, cache, inflow)
    own = SVector{3}(@view states[:, panel])
    residual_of = unknowns -> begin
        updated = ntuple(4) do side
            edge = cache.topology.cell_edges[panel, side]
            edge != 0 && cache.edge_donor[edge] == panel ? unknowns : donors[side]
        end
        cell_residual(unknowns, updated, panel, mesh, cache, shear)
    end
    current = residual_of(own)
    all(isfinite, current) || return false
    jacobian = ForwardDiff.jacobian(residual_of, own)
    step = try
        -(jacobian \ current)
    catch
        return false
    end
    all(isfinite, step) || return false

    largest = maximum(abs, step)
    largest > step_limit && (step = step * (step_limit / largest))
    damping = one(eltype(step))
    for _ in 1:16
        candidate = own + damping * step
        if _is_admissible(candidate, cache.closure)
            trial = residual_of(candidate)
            if all(isfinite, trial) && norm(trial) < norm(current)
                states[:, panel] .= candidate
                return true
            end
        end
        damping /= 2
    end
    return false
end

function _relax_panel!(states, panel, mesh, cache, inflow, shear)
    donors = _donor_states(states, panel, cache, inflow)
    guess = SVector{3}(@view states[:, panel])
    residual = (unknowns, _) -> cell_residual(unknowns, donors, panel, mesh, cache,
        shear)
    problem = NonlinearProblem(residual, guess, nothing)
    solution = NonlinearSolve.solve(problem, SimpleNewtonRaphson();
        abstol = 1e-12, reltol = 1e-12, maxiters = 30)
    all(isfinite, solution.u) || return false
    _is_admissible(solution.u, cache.closure) || return false
    # Accept on improvement rather than on the return code. The sweep only has
    # to hand Newton a better starting point, and discarding an unconverged but
    # improved iterate is worse than useless here: it freezes those cells at the
    # initial guess while their upstream neighbours advance, which leaves the
    # sweep with a larger residual than it started from and idempotent
    # thereafter.
    norm(residual(solution.u, nothing)) < norm(residual(guess, nothing)) || return false
    states[:, panel] .= solution.u
    return true
end

raw"""
    assemble_jacobian(states, mesh, cache)

Sparse Jacobian of the global residual.

The pattern is exactly the upwind cell-adjacency graph: because the velocity
gradients are frozen and the shear coefficient is lagged, panel `i`'s residual
depends only on its own state and on the states of the panels that donate flux
to it. Each edge therefore contributes one three-by-three block to each of its
two panels, differentiated once with respect to the donor state.
"""
function assemble_jacobian(states, mesh::Mesh, cache::SurfaceBoundaryLayerCache,
        inflow, shear; couple_shear::Bool = false)
    nfaces = mesh.nfaces
    rows = Int[]
    columns = Int[]
    values = eltype(states)[]
    blocks = 4 * nfaces + 2 * cache.topology.nedges
    sizehint!(rows, 9 * blocks)
    sizehint!(columns, 9 * blocks)
    sizehint!(values, 9 * blocks)

    function push_block!(block, panel, other)
        for local_row in 1:3, local_column in 1:3
            value = block[local_row, local_column]
            iszero(value) && continue
            push!(rows, 3 * (panel - 1) + local_row)
            push!(columns, 3 * (other - 1) + local_column)
            push!(values, value)
        end
    end

    for panel in 1:nfaces
        if !cache.active[panel]
            # Identity row: the residual above is constant for these panels.
            for component in 1:3
                push!(rows, 3 * (panel - 1) + component)
                push!(columns, 3 * (panel - 1) + component)
                push!(values, one(eltype(states)))
            end
            continue
        end
        donors = _donor_states(states, panel, cache, inflow)
        own = SVector{3}(@view states[:, panel])
        diffusive = _needs_neighbours(cache)
        neighbours = diffusive ? _neighbour_states(states, panel, cache) : nothing
        # Diagonal block: the source term, any edge whose donor is this panel,
        # and — when the diffusion term is on — its own side of every symmetric
        # edge flux.
        diagonal = ForwardDiff.jacobian(
            unknowns -> begin
                updated = ntuple(4) do side
                    edge = cache.topology.cell_edges[panel, side]
                    edge != 0 && cache.edge_donor[edge] == panel ? unknowns :
                    donors[side]
                end
                nearby = if diffusive
                    ntuple(4) do side
                        edge = cache.topology.cell_edges[panel, side]
                        edge == 0 && return neighbours[side]
                        other = cache.topology.edge_kind[edge] === :interior ?
                                edge_partner(cache.topology, edge, panel) : 0
                        (other == 0 || !cache.active[other]) ? unknowns :
                        neighbours[side]
                    end
                else
                    nothing
                end
                cell_residual(unknowns, updated, panel, mesh, cache,
                    _cell_shear(unknowns, shear, panel, mesh, cache, couple_shear),
                    nearby)
            end, own)
        push_block!(diagonal, panel, panel)

        # Every panel this one's residual reads: its donors, and — with
        # diffusion on — every active neighbour, since that term is symmetric.
        touched = Set{Int}()
        for side in 1:4
            edge = cache.topology.cell_edges[panel, side]
            edge == 0 && continue
            donor = cache.edge_donor[edge]
            donor != panel && donor != 0 && push!(touched, donor)
            diffusive && cache.topology.edge_kind[edge] === :interior || continue
            other = edge_partner(cache.topology, edge, panel)
            other != 0 && other != panel && cache.active[other] && push!(touched, other)
        end
        for other in touched
            block = ForwardDiff.jacobian(
                unknowns -> begin
                    updated = ntuple(4) do side
                        edge = cache.topology.cell_edges[panel, side]
                        edge != 0 && cache.edge_donor[edge] == other ? unknowns :
                        donors[side]
                    end
                    nearby = if diffusive
                        ntuple(4) do side
                            edge = cache.topology.cell_edges[panel, side]
                            edge == 0 && return neighbours[side]
                            partner = cache.topology.edge_kind[edge] === :interior ?
                                      edge_partner(cache.topology, edge, panel) : 0
                            partner == other ? unknowns : neighbours[side]
                        end
                    else
                        nothing
                    end
                    cell_residual(own, updated, panel, mesh, cache,
                        _cell_shear(own, shear, panel, mesh, cache, couple_shear),
                        nearby)
                end, SVector{3}(@view states[:, other]))
            push_block!(block, panel, other)
        end
    end
    return sparse(rows, columns, values, 3 * nfaces, 3 * nfaces)
end

raw"""
    global_residual(states, mesh, cache, inflow)

Residual of every panel, flattened to a vector of length ``3N``.
"""
function global_residual(states, mesh::Mesh, cache::SurfaceBoundaryLayerCache,
        inflow, shear; couple_shear::Bool = false)
    # The implicit-function path evaluates this with a primal state against a
    # cache built from dual inputs, which is the whole point of it: the
    # derivative enters through the edge velocity, not through the unknowns.
    element = promote_type(eltype(states), eltype(cache.speed), eltype(shear),
        eltype(mesh.areas))
    residual = Vector{element}(undef, 3 * mesh.nfaces)
    for panel in 1:mesh.nfaces
        own = SVector{3}(@view states[:, panel])
        # A panel outside the solved region keeps whatever state it was given
        # and contributes an identity row, so the system stays square and
        # nonsingular without that panel influencing anything.
        local_residual = if cache.active[panel]
            cell_residual(own, _donor_states(states, panel, cache, inflow), panel,
                mesh, cache, _cell_shear(own, shear, panel, mesh, cache, couple_shear),
                _needs_neighbours(cache) ? _neighbour_states(states, panel, cache) : nothing)
        else
            zero(SVector{3, element})
        end
        residual[(3 * (panel - 1) + 1):(3 * panel)] .= local_residual
    end
    return residual
end

raw"""
    _cell_shear(own, shear, panel, mesh, cache, couple)

The shear field a cell's residual should see.

With `couple == false` this is the lagged field itself, which is what the Newton
solve wants: the lag is deliberately outside the system so that the Jacobian
keeps the sparsity of the upwind graph.

With `couple == true` the panel's own entry is recomputed from its own state,
its donors' entries staying frozen. At a converged state this changes nothing —
the recomputation reproduces exactly what the outer loop converged to — so the
solution is untouched and only the *derivative* differs. It differs by about ten
per cent: freezing the lag entirely and differentiating the truncated iteration
bracket the true sensitivity from opposite sides, which is what says the lag is
a real part of the state and not a parameter. Coupling only the own-cell term
keeps the sparsity exactly as it was, and captures nearly all of it, because the
lag relaxes over far less than a panel — the donor's contribution enters through
`exp(-1.05 Δs/δ)`, of order a few per cent per cell.
"""
function _cell_shear(own, shear, panel, mesh, cache, couple::Bool)
    couple || return shear
    return _OverriddenShear(shear, panel,
        _panel_shear_of(own, shear, panel, mesh, cache))
end

raw"""
    solve_surface_boundary_layer(mesh, edge_velocity, kinematic_viscosity; kwargs...)

Solve the three-dimensional integral boundary layer over a panelled surface as
one globally coupled system.

This replaces the strip march of
[`solve_three_dimensional_boundary_layer`](@ref). The equations are the same —
two momentum-defect equations and a kinetic energy equation in a local Cartesian
surface basis — but they are discretised as a cell-centred finite volume over
the welded edge topology and solved everywhere at once:

```math
R_i=\frac{1}{A_i}\sum_{e\in\partial i}\ell_e\,\Phi_e^{(i)}
+g\left(U_i,\nabla q_i\right)=0 .
```

Nothing about the discretisation refers to mesh lines, so it does not matter
whether the external streamlines run along them. That is the defect the march
could not be repaired for.

The solve is a lagged-coefficient outer iteration wrapped around a damped sparse
Newton solve. The shear-stress coefficient is held fixed through each Newton
solve and refreshed between them; the Newton solve is then strictly monotone,
because the line search compares residuals of one and the same function. Letting
the lag move between Newton steps destroys that and the iteration stalls with
the shape factor drifting into separation.

The Jacobian sparsity is exactly the upwind adjacency graph, because the
velocity gradients are frozen and the shear coefficient lagged.

Nonlinear Gauss-Seidel sweeps in flow order (see [`flow_ordering`](@ref)) are
available as a globalization but are off by default: with the arclength estimate
of [`initial_states`](@ref) the Newton solve converges from the initial guess
directly, and a sweep that cannot converge a cell locally leaves it behind while
its upwind neighbours advance.

Accepts either a `Mesh` or a [`StructuredPanelGrid`](@ref); the strips are not
used.

# Keywords

- `outer_iterations`: maximum lagged-coefficient passes.
- `sweeps`: Gauss-Seidel passes before the first Newton solve.
- `newton_steps`: maximum damped Newton iterations per outer pass.
- `tolerance`: convergence threshold on the scaled residual norm.
- `shear_relaxation`: under-relaxation of the lagged shear coefficient.
- `step_limit`: largest move of any panel's transformed state in one Newton step.
- `gradient_bound`: cap on the velocity change one cell's gradients imply across
  it, relative to the local edge speed. See `_gradient_rescaling`.
- `implicit_derivative`: take derivatives through the implicit function theorem
  rather than by differentiating the iteration. The lagged shear coefficient is
  held at its converged value when the adjoint is formed, so the sensitivity is
  that of the frozen-coefficient system; the lag relaxes over far less than a
  panel, which is why it was lagged in the first place. Set `false` to
  differentiate the iteration itself, which is slower and tolerance-dependent
  but exact for the lag.
- `topology`, `metrics`, `cache`: supply precomputed geometry to avoid rebuilding
  it across a coupling loop or a derivative sweep.
"""
function solve_surface_boundary_layer(mesh::Mesh, edge_velocity::AbstractMatrix,
        kinematic_viscosity::Real; rho::Real = SETTINGS.rho, reference = (0, 0, 0),
        closure::ThreeDimensionalClosure = ThreeDimensionalClosure(),
        topology = nothing, metrics = nothing, cache = nothing,
        minimum_edge_speed = nothing, weld_tolerance = nothing, gradient_bound = 1.0,
        sweeps::Integer = 0,
        outer_iterations::Integer = 25, newton_steps::Integer = 60,
        tolerance::Real = 1e-10, shear_relaxation::Real = 0.5, step_limit::Real = 1.0,
        implicit_derivative::Bool = true, active = nothing, diffusion::Real = 0.0,
        upwind_smoothing::Real = 0.0, initial = nothing)
    rho > 0 || throw(ArgumentError("rho must be positive"))
    sweeps >= 0 || throw(ArgumentError("sweeps must be nonnegative"))
    newton_steps >= 0 || throw(ArgumentError("newton_steps must be nonnegative"))
    outer_iterations >= 1 ||
        throw(ArgumentError("outer_iterations must be at least one"))
    0 < shear_relaxation <= 1 ||
        throw(ArgumentError("shear_relaxation must lie in (0, 1]"))
    step_limit > 0 || throw(ArgumentError("step_limit must be positive"))

    settings = (; closure, topology, metrics, minimum_edge_speed, weld_tolerance,
        gradient_bound, sweeps, outer_iterations, newton_steps, tolerance,
        shear_relaxation, step_limit, initial)

    if !implicit_derivative
        problem = isnothing(cache) ?
                  build_surface_cache(mesh, edge_velocity, kinematic_viscosity;
            closure, topology, metrics, minimum_edge_speed, weld_tolerance,
            gradient_bound) : cache
        solution = _converge_states(mesh, problem, kinematic_viscosity, settings)
        return assemble_surface_result(mesh, problem, solution, rho, reference)
    end

    # The converged state satisfies R(y, x) = 0, so its derivative comes from
    # the implicit function theorem rather than from differentiating the
    # iteration. That makes the sensitivity independent of how many Newton
    # steps were taken and of the tolerance, and costs one back-substitution per
    # seed against one full solve per seed. The Jacobian the adjoint needs is
    # the one Newton already assembles, so the two are the same piece of work.
    lagged = Ref{Any}(nothing)
    diagnostics = Ref{Any}(nothing)
    shape = size(edge_velocity)

    function converge(inputs, _)
        velocity = reshape(inputs, shape)
        problem = build_surface_cache(mesh, velocity, kinematic_viscosity; closure,
            topology, metrics, minimum_edge_speed, weld_tolerance, gradient_bound,
        active, diffusion, upwind_smoothing)
        solution = _converge_states(mesh, problem, kinematic_viscosity, settings)
        # Called once, with the primal inputs. Keeping the lagged shear and the
        # diagnostics from that call is what lets the residual below be a
        # function of the state alone.
        lagged[] = solution.shear_coefficient
        diagnostics[] = solution
        return vec(solution.states)
    end

    function residual_of(unknowns, inputs, _)
        velocity = reshape(inputs, shape)
        problem = build_surface_cache(mesh, velocity, kinematic_viscosity; closure,
            topology, metrics, minimum_edge_speed, weld_tolerance, gradient_bound,
        active, diffusion, upwind_smoothing)
        return global_residual(reshape(unknowns, 3, mesh.nfaces), mesh, problem,
            inflow_states(mesh, problem), lagged[]; couple_shear = true)
    end

    function jacobian_of(_, unknowns, inputs, _)
        velocity = reshape(inputs, shape)
        problem = build_surface_cache(mesh, velocity, kinematic_viscosity; closure,
            topology, metrics, minimum_edge_speed, weld_tolerance, gradient_bound,
        active, diffusion, upwind_smoothing)
        return assemble_jacobian(reshape(unknowns, 3, mesh.nfaces), mesh, problem,
            inflow_states(mesh, problem), lagged[]; couple_shear = true)
    end

    # Factorise once and back-substitute every seed. ImplicitAD hands this a
    # MATRIX right-hand side, one column per seed, which is exactly what `\`
    # after `lu` is for; a LinearSolve cache buys nothing here because the
    # adjoint is solved once rather than repeatedly.
    converged = ImplicitAD.implicit(converge, residual_of, vec(edge_velocity);
        drdy = jacobian_of, lsolve = (matrix, right) -> lu(matrix) \ right)

    problem = build_surface_cache(mesh, edge_velocity, kinematic_viscosity; closure,
        topology, metrics, minimum_edge_speed, weld_tolerance, gradient_bound,
        active, diffusion, upwind_smoothing)
    reported = diagnostics[]
    # Promoting the frozen lag to the state's element type gives it zero
    # partials, which is precisely the frozen-coefficient sensitivity described
    # above rather than an accident of the container being homogeneous.
    element = eltype(converged)
    solution = SurfaceBoundaryLayerSolution(reshape(converged, 3, mesh.nfaces),
        convert(Vector{element}, lagged[]), convert(element, reported.residual_norm),
        reported.sweeps, reported.newton_steps, reported.converged)
    return assemble_surface_result(mesh, problem, solution, rho, reference)
end

# The iteration itself, on whatever number type it is handed. Split out so that
# the differentiated interface above can run it once on primal values.
function _converge_states(mesh::Mesh, problem::SurfaceBoundaryLayerCache,
        kinematic_viscosity, settings)
    order = flow_ordering(problem)
    inflow = inflow_states(mesh, problem)
    # A warm start may arrive as plain numbers while this solve carries dual
    # numbers, which is exactly what happens when a coupling loop reuses the
    # base state to seed a differentiated one, so promote rather than copy.
    states = if isnothing(settings.initial)
        initial_states(mesh, problem, order, inflow)
    else
        size(settings.initial) == (3, mesh.nfaces) || throw(DimensionMismatch(
            "initial must have size (3, mesh.nfaces)"))
        convert(Matrix{promote_type(eltype(settings.initial), eltype(problem.speed))},
            settings.initial)
    end

    scale = _residual_scale(mesh, problem)
    weights = _residual_weights(mesh)
    shear = _lagged_shear(states, mesh, problem, order)
    completed_sweeps = 0
    for _ in 1:settings.sweeps
        for panel in order
            _march_panel!(states, panel, mesh, problem, inflow, shear,
                settings.step_limit)
            shear[panel] = _panel_shear(states, shear, panel, mesh, problem)
        end
        completed_sweeps += 1
        _weighted_norm(global_residual(states, mesh, problem, inflow, shear),
            weights) / scale < settings.tolerance && break
    end

    completed_newton = 0
    residual = global_residual(states, mesh, problem, inflow, shear)
    for _ in 1:settings.outer_iterations
        completed_newton += _newton_solve!(states, residual, mesh, problem, inflow,
            shear, settings.closure, scale, settings.tolerance, settings.newton_steps,
            settings.step_limit, weights)
        refreshed = _lagged_shear(states, mesh, problem, order)
        drift = maximum(abs, refreshed .- shear) /
                max(maximum(abs, refreshed), eps(Float64))
        @. shear = shear + settings.shear_relaxation * (refreshed - shear)
        residual = global_residual(states, mesh, problem, inflow, shear)
        _weighted_norm(residual, weights) / scale < settings.tolerance &&
            drift < settings.tolerance && break
    end

    residual_norm = _weighted_norm(residual, weights) / scale
    return SurfaceBoundaryLayerSolution(states, shear, residual_norm,
        completed_sweeps, completed_newton, residual_norm < settings.tolerance)
end

raw"""
Damped Newton on the flux residual with the shear coefficient held fixed, so
that the line search compares like with like and the iteration is monotone.
`residual` is updated in place of the caller's value through the return of the
step count; the caller re-evaluates once the lag moves.

# Why this is not `NonlinearSolve`

`NonlinearSolve` *is* used, for the three-by-three cell problems in
[`_relax_panel!`](@ref). It is not used for the global system, and the reason is
measured rather than assumed. Run against this same residual on the KVLCC2
forebody, its solvers give:

| solver | scaled residual | state admissible |
|---|---|---|
| this function | 2.691 | yes |
| `NewtonRaphson` | 9.078 | no |
| `TrustRegion` | 9.078 | no |
| `RobustMultiNewton` | 1.983 | **no** |
| `Broyden` | 9.078 | no |

`RobustMultiNewton` reaches a *lower* residual than this function and the result
is still useless, because it leaves the region where the closure relations are
defined — the shape factor runs past the bound `_is_admissible` enforces. A root
of the residual outside the closure's validity is not a solution of the boundary
layer.

That constraint is not expressible through the `NonlinearSolve` interface: it is
a bound on the *state*, enforced here by shortening each panel's step
individually, which is a reparameterisation of the step rather than a scalar
line search along a fixed direction.

The clean fix is not a different solver but a different parameterisation. Lokatt
and Eller bound the kinematic shape factor by construction, taking
``H_k=1+k_1/(1+(k_2u_1)^2)``, so that no iterate can leave validity and any
off-the-shelf solver applies. Adopting that would make `NonlinearSolve` usable
here directly, at the cost of changing the state variables everything else is
written against.
"""
function _newton_solve!(states, residual, mesh, cache, inflow, shear, closure, scale,
        tolerance, newton_steps, step_limit, weights)
    taken = 0
    current = residual
    linear = Ref{Any}(nothing)
    for _ in 1:newton_steps
        _weighted_norm(current, weights) / scale < tolerance && break
        jacobian = assemble_jacobian(states, mesh, cache, inflow, shear)
        step = _newton_step(linear, jacobian, current)
        isnothing(step) && break
        accepted = false
        damping = one(eltype(states))
        increment = reshape(step, 3, mesh.nfaces)
        candidate = similar(states)
        for _ in 1:16
            # Try the damped Newton step as it stands before capping how far any
            # one panel may move. The cap exists so that a few panels asking for
            # enormous steps cannot wreck the direction for everyone, but it is
            # not free: it changes the direction, and once the iterate is close
            # the uncapped step is the better one. Capping unconditionally left
            # a measured third of the available residual reduction unused and
            # stalled the solve; releasing it unconditionally sends the early
            # iterations somewhere worse. Trying both and keeping whichever
            # descends costs one extra residual evaluation per damping level.
            for limit in (typemax(step_limit), step_limit)
                _clipped_candidate!(candidate, states, increment, damping, closure,
                    limit)
                trial = global_residual(candidate, mesh, cache, inflow, shear)
                if all(isfinite, trial) &&
                   _weighted_norm(trial, weights) < _weighted_norm(current, weights)
                    states .= candidate
                    current = trial
                    accepted = true
                    break
                end
            end
            accepted && break
            damping /= 2
        end
        taken += 1
        accepted || break
    end
    return taken
end

raw"""
    _clipped_candidate!(candidate, states, increment, damping, closure, step_limit)

Apply `damping * increment` to every panel, shortening it *per panel* so that
none leaves the closure's validity region and none moves further than
`step_limit` in the transformed variables.

Clipping per panel rather than rejecting the whole step is what makes the solve
usable on a hull. The stem carries a handful of tiny, strongly skewed panels
sitting on the stagnation line, and their Newton steps are enormous next to the
rest of the surface; an all-or-nothing admissibility test lets any one of them
veto the step for the other several hundred, so the line search collapses to a
damping of order 1e-5 and the iteration never moves at all.

The line search still requires the full residual to fall, so the clipping only
changes the path taken, never the solution accepted. Near convergence every
panel's step is small, no clipping is active, and the step is exactly Newton's.
"""
function _clipped_candidate!(candidate, states, increment, damping, closure,
        step_limit)
    for panel in axes(states, 2)
        own = SVector{3}(@view states[:, panel])
        move = damping * SVector{3}(@view increment[:, panel])
        largest = maximum(abs, move)
        largest > step_limit && (move = move * (step_limit / largest))
        # Shorten until the panel is admissible. The transformed variables keep
        # the iterate inside the closure for all but the largest excursions, so
        # this loop almost always exits immediately.
        for _ in 1:20
            _is_admissible(own + move, closure) && break
            move = move / 2
        end
        candidate[:, panel] .= own .+ move
    end
    return candidate
end

solve_surface_boundary_layer(grid::StructuredPanelGrid, args...; kwargs...) =
    solve_surface_boundary_layer(grid.mesh, args...; kwargs...)

# Residuals carry the units of the source terms, whose dominant scale is the
# skin-friction coefficient. Normalising by it makes the tolerance mesh- and
# speed-independent.
function _residual_scale(mesh, cache)
    total = zero(eltype(cache.speed))
    for panel in 1:mesh.nfaces
        total += cache.speed[panel]
    end
    return max(sqrt(mesh.nfaces) * 1e-3, sqrt(eps(Float64)) * total)
end

raw"""
    _residual_weights(mesh)

Per-panel weights making the residual norm a discrete surface integral rather
than a raw sum over cells.

A cell residual is a flux divergence, so it scales inversely with the cell size;
weighting by ``\sqrt{A_i/\bar A}`` recovers the ``L^2`` norm of the underlying
field. Without it the smallest panel dictates the convergence test. On the
KVLCC2 hull the areas span a factor of 150 and a handful of slivers at the stem
hold nine tenths of the unweighted squared residual, so the line search spends
the whole solve on cells covering a thousandth of the surface.

Row scaling by a fixed diagonal leaves the Newton direction untouched — the
Jacobian and the residual scale together — so this changes only the merit
function the line search descends and the quantity the tolerance is applied to.
"""
function _residual_weights(mesh::Mesh)
    mean_area = sum(mesh.areas) / mesh.nfaces
    return sqrt.(mesh.areas ./ mean_area)
end

# Norm of a stacked cell residual under the area weighting above.
function _weighted_norm(residual, weights)
    total = zero(eltype(residual))
    for panel in eachindex(weights)
        base = 3 * (panel - 1)
        for component in 1:3
            total += (weights[panel] * residual[base + component])^2
        end
    end
    return sqrt(total)
end

# The shear-stress coefficient is advanced along the flow with the exact
# relaxation, using the dominant inflow edge to set the step. It is lagged out
# of the Newton system on purpose: it relaxes over a length far shorter than a
# panel, so including it would stiffen the Jacobian for no accuracy gain.
#
# One panel at a time, so that a Gauss-Seidel sweep can refresh the lag as it
# goes. Solving a cell against one shear field and then measuring its residual
# against another is not a fixed point of anything, and it costs the sweep its
# march-like behaviour.
raw"""
    _OverriddenShear(base, index, value)

The lagged shear field with one entry replaced.

Used to give a panel's own shear coefficient a derivative with respect to its
own state while leaving its donors' frozen. Lazy because the alternative —
copying the whole field once per panel — is quadratic in the panel count.
"""
struct _OverriddenShear{V, T} <: AbstractVector{T}
    base::V
    index::Int
    value::T
end

Base.size(shear::_OverriddenShear) = size(shear.base)
Base.@propagate_inbounds function Base.getindex(shear::_OverriddenShear, index::Int)
    return index == shear.index ? shear.value :
           convert(typeof(shear.value), shear.base[index])
end

function _panel_shear(states, shear, panel, mesh, cache)
    return _panel_shear_of(SVector{3}(@view states[:, panel]), shear, panel, mesh,
        cache)
end

function _panel_shear_of(own, shear, panel, mesh, cache)
    element_type = eltype(own)
    state = panel_fluxes(own, cache.speed[panel], cache.cosine[panel],
        cache.sine[panel], zero(element_type), cache.viscosity, cache.closure)
    upstream = zero(element_type)
    step = sqrt(mesh.areas[panel])
    strongest = zero(element_type)
    for side in 1:4
        edge = cache.topology.cell_edges[panel, side]
        edge == 0 && continue
        donor = cache.edge_donor[edge]
        (donor == 0 || donor == panel) && continue
        strength = _inflow_strength(cache, panel, edge)
        if strength > strongest
            area_vector = edge_area_vector_for(cache.topology, cache.metrics, edge,
                panel)
            strongest = strength
            upstream = shear[donor]
            step = mesh.areas[panel] /
                   max(hypot(area_vector[1], area_vector[2]), eps(Float64))
        end
    end
    equilibrium = equilibrium_shear_coefficient(state.shape_factor, state.kinetic,
        state.slip)
    return upstream > 0 ? advance_shear_coefficient(upstream, state, step) : equilibrium
end

function _lagged_shear(states, mesh, cache, order)
    shear = zeros(eltype(states), mesh.nfaces)
    for panel in order
        shear[panel] = _panel_shear(states, shear, panel, mesh, cache)
    end
    return shear
end

raw"""
    assemble_surface_result(mesh, cache, solution, rho, reference)

Package a converged surface solve as a
[`ThreeDimensionalBoundaryLayerResult`](@ref).

The same per-panel result type as the strip march, deliberately: nothing in it
is march-specific, and `viscous_maneuvering_correction` consumes only
`transpiration_velocity`, `force[2]`, `moment[3]` and `separated`. Reusing it
makes this solver a drop-in.

The transpiration is the surface divergence of the mass defect,

```math
v_n=\frac{1}{A_i}\sum_{e\in\partial i}
\ell_e\left(q_e\boldsymbol\delta^*\right)\cdot\boldsymbol m_e ,
```

taken over the *full* displacement vector. Its crossflow component is the
physical content the strip march discarded entirely.
"""
function assemble_surface_result(mesh::Mesh, cache::SurfaceBoundaryLayerCache,
        solution::SurfaceBoundaryLayerSolution, rho, reference)
    element_type = eltype(solution.states)
    nfaces = mesh.nfaces
    edge_speed = zeros(element_type, nfaces)
    momentum_thickness = zeros(element_type, nfaces)
    shape_factor = zeros(element_type, nfaces)
    kinetic_shape_factor = zeros(element_type, nfaces)
    crossflow_angle = zeros(element_type, nfaces)
    displacement = zeros(element_type, nfaces)
    crossflow_displacement = zeros(element_type, nfaces)
    friction = zeros(element_type, nfaces)
    wall_shear = zeros(element_type, nfaces, 3)
    transpiration = zeros(element_type, nfaces)
    separated = falses(nfaces)

    states = Vector{Any}(undef, nfaces)
    for panel in 1:nfaces
        own = SVector{3}(@view solution.states[:, panel])
        state = panel_fluxes(own, cache.speed[panel], cache.cosine[panel],
            cache.sine[panel], solution.shear_coefficient[panel], cache.viscosity,
            cache.closure)
        states[panel] = state
        edge_speed[panel] = cache.speed[panel]
        momentum_thickness[panel] = state.momentum_thickness
        shape_factor[panel] = state.shape_factor
        kinetic_shape_factor[panel] = state.kinetic
        crossflow_angle[panel] = state.crossflow_angle
        displacement[panel] = state.displacement[1]
        crossflow_displacement[panel] = state.displacement[2]
        friction[panel] = state.friction_magnitude
        separated[panel] = state.shape_factor >= cache.closure.separation_shape_factor

        scale = separated[panel] ? zero(element_type) :
                0.5 * rho * cache.speed[panel]^2
        for axis in 1:3
            wall_shear[panel, axis] = scale *
                                      (state.friction[1] *
                                       cache.metrics.tangent[panel, axis] +
                                       state.friction[2] *
                                       cache.metrics.binormal[panel, axis])
        end
    end

    # Surface divergence of the mass defect, upwinded consistently with the
    # conservation equations.
    for panel in 1:nfaces
        total = zero(element_type)
        for side in 1:4
            edge = cache.topology.cell_edges[panel, side]
            edge == 0 && continue
            donor = cache.edge_donor[edge]
            source = donor == 0 ? panel : donor
            state = states[source]
            defect = @SVector [cache.speed[source] * state.displacement[1],
                cache.speed[source] * state.displacement[2]]
            if source != panel
                defect = rotation_into(cache.topology, cache.metrics, edge, panel) *
                         defect
            end
            area_vector = edge_area_vector_for(cache.topology, cache.metrics, edge,
                panel)
            total += defect[1] * area_vector[1] + defect[2] * area_vector[2]
        end
        transpiration[panel] = separated[panel] ? zero(element_type) :
                               total / mesh.areas[panel]
    end

    force = zeros(element_type, 3)
    moment = zeros(element_type, 3)
    for panel in 1:nfaces
        (separated[panel] || !cache.active[panel]) && continue
        panel_force = @view(wall_shear[panel, :]) .* mesh.areas[panel]
        force .+= panel_force
        moment .+= cross(@view(mesh.centers[panel, :]) .- reference, panel_force)
    end

    return ThreeDimensionalBoundaryLayerResult(edge_speed, momentum_thickness,
        shape_factor, kinetic_shape_factor, crossflow_angle,
        solution.shear_coefficient, displacement, crossflow_displacement, friction,
        wall_shear, transpiration, separated, .!separated, force, moment,
        (converged = solution.converged, residual_norm = solution.residual_norm,
            sweeps = solution.sweeps, newton_steps = solution.newton_steps,
            states = solution.states))
end
