import NonlinearSolve

using LinearAlgebra: cross, dot, lu, norm
using NonlinearSolve: NonlinearProblem, SimpleNewtonRaphson
using SparseArrays: SparseMatrixCSC, sparse

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
        weld_tolerance = nothing, gradient_bound::Real = 1.0)
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
        along = dot(velocity, @view(geometry.tangent[panel, :]))
        across = dot(velocity, @view(geometry.binormal[panel, :]))
        magnitude = max(hypot(along, across), floor_speed)
        speed[panel] = magnitude
        cosine[panel] = along / magnitude
        sine[panel] = across / magnitude
    end

    gradients = _least_squares_gradients(mesh, surface, geometry, speed, cosine, sine,
        element_type, gradient_bound)
    donor, influx = _upwind_donors(mesh, surface, geometry, speed, cosine, sine)

    return SurfaceBoundaryLayerCache(surface, geometry, speed, cosine, sine,
        gradients, donor, influx, convert(element_type, kinematic_viscosity), closure)
end

# Weighted linear least squares for the in-surface gradient of the edge velocity
# and its magnitude. The moment matrix depends only on geometry, so the whole
# stencil could be cached; it is rebuilt here because the metrics carry
# derivatives and the cost is negligible next to the Newton solve.
function _least_squares_gradients(mesh, topology, metrics, speed, cosine, sine,
        element_type, bound)
    gradient_type = NamedTuple{
        (:speed_x, :speed_y, :u_x, :u_y, :v_x, :v_y), NTuple{6, element_type}}
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
            offset3 = @view(mesh.centers[other, :]) .- @view(mesh.centers[panel, :])
            offset = @SVector [dot(offset3, @view(metrics.tangent[panel, :])),
                dot(offset3, @view(metrics.binormal[panel, :]))]
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
    largest = max(abs(speed_gradient[1]), abs(speed_gradient[2]), abs(u_gradient[1]),
        abs(u_gradient[2]), abs(v_gradient[1]), abs(v_gradient[2]))
    scaled = largest * size / max(speed, eps(typeof(speed)))
    return scaled > bound ? bound / scaled : one(scaled)
end

function _upwind_donors(mesh, topology, metrics, speed, cosine, sine)
    donor = zeros(Int, topology.nedges)
    influx = zeros(eltype(speed), topology.nedges)
    for edge in 1:topology.nedges
        left = topology.edge_cells[edge, 1]
        right = topology.edge_cells[edge, 2]
        normal = (metrics.edge_area_vector[edge, 1, 1],
            metrics.edge_area_vector[edge, 1, 2])
        if right == 0
            # On a boundary edge, outflow extrapolates from the owning cell and
            # inflow is marked with a zero donor so that a prescribed state is
            # supplied instead. Without that the flux divergence of a leading
            # cell would cancel by the geometric conservation law and the
            # momentum equation would have no solution.
            velocity = (speed[left] * cosine[left], speed[left] * sine[left])
            outward = velocity[1] * normal[1] + velocity[2] * normal[2]
            donor[edge] = outward >= 0 ? left : 0
            influx[edge] = abs(outward)
            continue
        end
        left_velocity = (speed[left] * cosine[left], speed[left] * sine[left])
        rotation = rotation_into(topology, metrics, edge, left)
        right_velocity = rotation *
                         @SVector([speed[right] * cosine[right],
            speed[right] * sine[right]])
        left_weight = one(eltype(speed)) /
                      max(norm(@view(metrics.edge_midpoint[edge, :]) .-
                               @view(mesh.centers[left, :])), eps(Float64))
        right_weight = one(eltype(speed)) /
                       max(norm(@view(metrics.edge_midpoint[edge, :]) .-
                                @view(mesh.centers[right, :])), eps(Float64))
        blended = ((left_weight * left_velocity[1] + right_weight * right_velocity[1]),
            (left_weight * left_velocity[2] + right_weight * right_velocity[2])) ./
                  (left_weight + right_weight)
        outward = blended[1] * normal[1] + blended[2] * normal[2]
        donor[edge] = outward >= 0 ? left : right
        influx[edge] = abs(outward)
    end
    return donor, influx
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
        cache::SurfaceBoundaryLayerCache, shear)
    donor = cache.edge_donor[edge]
    # An inflow boundary edge carries the prescribed state, evaluated in the
    # receiving panel's own frame.
    source_panel = donor == 0 ? panel : donor
    state = panel_fluxes(donor_state, cache.speed[source_panel],
        cache.cosine[source_panel], cache.sine[source_panel],
        shear[source_panel], cache.viscosity, cache.closure)
    momentum = state.momentum
    energy = state.energy
    if donor != panel && donor != 0
        rotation = rotation_into(cache.topology, cache.metrics, edge, panel)
        momentum = rotation * momentum * transpose(rotation)
        energy = rotation * energy
    end
    area_vector = edge_area_vector_for(cache.topology, cache.metrics, edge, panel)
    normal = @SVector [area_vector[1], area_vector[2]]
    transported = momentum * normal
    return @SVector [transported[1], transported[2], dot(energy, normal)]
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
        cache::SurfaceBoundaryLayerCache, shear)
    flux_sum = zero(SVector{3, eltype(own_state)})
    for side in 1:4
        edge = cache.topology.cell_edges[panel, side]
        edge == 0 && continue
        flux_sum = flux_sum +
                   edge_flux_contribution(donor_states[side], edge, panel, cache,
            shear)
    end
    state = panel_fluxes(own_state, cache.speed[panel], cache.cosine[panel],
        cache.sine[panel], shear[panel], cache.viscosity, cache.closure)
    sources = station_sources(state, (speed = cache.speed[panel],), cache.gradients[panel])
    return flux_sum ./ mesh.areas[panel] .+ sources
end

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
            offset = @view(mesh.centers[panel, :]) .- @view(mesh.centers[donor, :])
            step = abs(dot(offset, @view(cache.metrics.tangent[panel, :])) *
                       cache.cosine[panel] +
                       dot(offset, @view(cache.metrics.binormal[panel, :])) *
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
        receiver = edge_partner(topology, edge, donor)
        # A side the flow runs parallel to transports nothing, so it must not
        # constrain the order: on a plate meshed square to the stream every
        # transverse side would otherwise chain the rows together and make the
        # sweep march across the flow as well as along it.
        _inflow_strength(cache, receiver, edge) > 0 || continue
        push!(downstream[donor], receiver)
        upstream_count[receiver] += 1
    end
    ready = [panel for panel in 1:nfaces if upstream_count[panel] == 0]
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
        append!(order, (panel for panel in 1:nfaces if !visited[panel]))
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
        inflow, shear)
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
        donors = _donor_states(states, panel, cache, inflow)
        own = SVector{3}(@view states[:, panel])
        # Diagonal block: the source term plus any edge whose donor is this
        # panel.
        diagonal = ForwardDiff.jacobian(
            unknowns -> begin
                updated = ntuple(4) do side
                    edge = cache.topology.cell_edges[panel, side]
                    edge != 0 && cache.edge_donor[edge] == panel ? unknowns :
                    donors[side]
                end
                cell_residual(unknowns, updated, panel, mesh, cache, shear)
            end, own)
        push_block!(diagonal, panel, panel)

        for side in 1:4
            edge = cache.topology.cell_edges[panel, side]
            edge == 0 && continue
            donor = cache.edge_donor[edge]
            (donor == panel || donor == 0) && continue
            block = ForwardDiff.jacobian(
                unknowns -> begin
                    updated = ntuple(4) do other_side
                        other_edge = cache.topology.cell_edges[panel, other_side]
                        other_edge == edge ? unknowns : donors[other_side]
                    end
                    cell_residual(own, updated, panel, mesh, cache, shear)
                end, SVector{3}(@view states[:, donor]))
            push_block!(block, panel, donor)
        end
    end
    return sparse(rows, columns, values, 3 * nfaces, 3 * nfaces)
end

raw"""
    global_residual(states, mesh, cache, inflow)

Residual of every panel, flattened to a vector of length ``3N``.
"""
function global_residual(states, mesh::Mesh, cache::SurfaceBoundaryLayerCache,
        inflow, shear)
    residual = Vector{eltype(states)}(undef, 3 * mesh.nfaces)
    for panel in 1:mesh.nfaces
        own = SVector{3}(@view states[:, panel])
        local_residual = cell_residual(own,
            _donor_states(states, panel, cache, inflow), panel, mesh, cache, shear)
        residual[(3 * (panel - 1) + 1):(3 * panel)] .= local_residual
    end
    return residual
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
        initial = nothing)
    rho > 0 || throw(ArgumentError("rho must be positive"))
    sweeps >= 0 || throw(ArgumentError("sweeps must be nonnegative"))
    newton_steps >= 0 || throw(ArgumentError("newton_steps must be nonnegative"))
    outer_iterations >= 1 ||
        throw(ArgumentError("outer_iterations must be at least one"))
    0 < shear_relaxation <= 1 ||
        throw(ArgumentError("shear_relaxation must lie in (0, 1]"))
    step_limit > 0 || throw(ArgumentError("step_limit must be positive"))

    problem = isnothing(cache) ?
              build_surface_cache(mesh, edge_velocity, kinematic_viscosity; closure,
        topology, metrics, minimum_edge_speed, weld_tolerance, gradient_bound) : cache
    order = flow_ordering(problem)
    inflow = inflow_states(mesh, problem)
    # A warm start may arrive as plain numbers while this solve carries dual
    # numbers, which is exactly what happens when a coupling loop reuses the
    # base state to seed a differentiated one, so promote rather than copy.
    states = if isnothing(initial)
        initial_states(mesh, problem, order, inflow)
    else
        size(initial) == (3, mesh.nfaces) || throw(DimensionMismatch(
            "initial must have size (3, mesh.nfaces)"))
        convert(Matrix{promote_type(eltype(initial), eltype(problem.speed))}, initial)
    end

    scale = _residual_scale(mesh, problem)
    weights = _residual_weights(mesh)
    shear = _lagged_shear(states, mesh, problem, order)
    completed_sweeps = 0
    for _ in 1:sweeps
        for panel in order
            _relax_panel!(states, panel, mesh, problem, inflow, shear)
            shear[panel] = _panel_shear(states, shear, panel, mesh, problem)
        end
        completed_sweeps += 1
        _weighted_norm(global_residual(states, mesh, problem, inflow, shear),
            weights) / scale < tolerance && break
    end

    completed_newton = 0
    residual = global_residual(states, mesh, problem, inflow, shear)
    for _ in 1:outer_iterations
        completed_newton += _newton_solve!(states, residual, mesh, problem, inflow,
            shear, closure, scale, tolerance, newton_steps, step_limit, weights)
        refreshed = _lagged_shear(states, mesh, problem, order)
        drift = maximum(abs, refreshed .- shear) /
                max(maximum(abs, refreshed), eps(Float64))
        @. shear = shear + shear_relaxation * (refreshed - shear)
        residual = global_residual(states, mesh, problem, inflow, shear)
        _weighted_norm(residual, weights) / scale < tolerance && drift < tolerance &&
            break
    end

    residual_norm = _weighted_norm(residual, weights) / scale
    solution = SurfaceBoundaryLayerSolution(states, shear, residual_norm,
        completed_sweeps, completed_newton, residual_norm < tolerance)
    return assemble_surface_result(mesh, problem, solution, rho, reference)
end

# Damped Newton on the flux residual with the shear coefficient held fixed, so
# that the line search compares like with like and the iteration is monotone.
# `residual` is updated in place of the caller's value through the return of the
# step count; the caller re-evaluates once the lag moves.
function _newton_solve!(states, residual, mesh, cache, inflow, shear, closure, scale,
        tolerance, newton_steps, step_limit, weights)
    taken = 0
    current = residual
    for _ in 1:newton_steps
        _weighted_norm(current, weights) / scale < tolerance && break
        jacobian = assemble_jacobian(states, mesh, cache, inflow, shear)
        step = try
            -(lu(jacobian) \ current)
        catch
            break
        end
        all(isfinite, step) || break
        accepted = false
        damping = one(eltype(states))
        increment = reshape(step, 3, mesh.nfaces)
        candidate = similar(states)
        for _ in 1:16
            _clipped_candidate!(candidate, states, increment, damping, closure,
                step_limit)
            trial = global_residual(candidate, mesh, cache, inflow, shear)
            if all(isfinite, trial) &&
               _weighted_norm(trial, weights) < _weighted_norm(current, weights)
                states .= candidate
                current = trial
                accepted = true
                break
            end
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
function _panel_shear(states, shear, panel, mesh, cache)
    element_type = eltype(states)
    own = SVector{3}(@view states[:, panel])
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
        separated[panel] && continue
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
