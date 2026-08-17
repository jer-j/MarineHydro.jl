raw"""
    SimultaneousCouplingResult

Converged state of a strongly coupled viscous-inviscid solve.

`layer` is the boundary-layer solution, `sources` the transpiration source
strengths, `edge_velocity` the converged external velocity the layer actually
saw, and `residual_norm` the scaled norm of the joint system at convergence.
"""
struct SimultaneousCouplingResult{L, T, M}
    layer::L
    sources::Vector{T}
    edge_velocity::M
    residual_norm::T
    newton_steps::Int
    converged::Bool
end

raw"""
    indirect_gradient_operator(green_functions, mesh, wavenumber = 0)

Dense operator ``G`` with ``\nabla\phi = G\sigma``, one row block per panel.

[`evaluate_indirect_potential_gradient`](@ref) is linear in the sources, so the
whole map can be assembled once and reused. That is what makes a strongly
coupled solve affordable: the inviscid influence enters the joint Jacobian as a
constant matrix rather than as a solve to be repeated per Newton step.

Returned as `(nfaces, 3, nfaces)`, indexed by target panel, axis, source panel.
"""
function indirect_gradient_operator(green_functions, mesh::Mesh,
        wavenumber::Real = 0)
    element_type = eltype(mesh.centers)
    operator = zeros(element_type, mesh.nfaces, 3, mesh.nfaces)
    free_surface = zero(mesh.centers[1, 3])
    for target_index in 1:mesh.nfaces
        target = element(mesh, target_index)
        one_coordinate = one(mesh.centers[target_index, 3])
        jump = abs(mesh.centers[target_index, 3] - free_surface) < 1e-8 ?
               one_coordinate : one_coordinate / 2
        target_normal = normal(target)
        for source_index in 1:mesh.nfaces
            influence = integral_gradient(green_functions, target,
                element(mesh, source_index), wavenumber;
                with_respect_to_first_variable = true)
            for axis in 1:3
                operator[target_index, axis, source_index] = -real(influence[axis]) /
                                                             (2τ̅)
            end
        end
        for axis in 1:3
            operator[target_index, axis, target_index] += jump * target_normal[axis]
        end
    end
    return operator
end

# Tangentially projected correction velocity from a source distribution.
function _correction_velocity(operator, mesh::Mesh, sources)
    element_type = promote_type(eltype(operator), eltype(sources))
    velocity = zeros(element_type, mesh.nfaces, 3)
    for panel in 1:mesh.nfaces
        normal_component = zero(element_type)
        for axis in 1:3
            total = zero(element_type)
            for source in eachindex(sources)
                total += operator[panel, axis, source] * sources[source]
            end
            velocity[panel, axis] = total
            normal_component += total * mesh.normals[panel, axis]
        end
        for axis in 1:3
            velocity[panel, axis] -= normal_component * mesh.normals[panel, axis]
        end
    end
    return velocity
end

raw"""
    solve_simultaneous_coupling(
        grid, rigid_edge_velocity, kinematic_viscosity;
        green_functions = (Rankine(), RankineReflected()),
        rho = SETTINGS.rho, closure = ThreeDimensionalClosure(), kwargs...
    )

Solve the boundary layer and the transpiration potential as **one** Newton
system.

The weak alternative already in [`viscous_maneuvering_correction`](@ref)
alternates: converge the layer against a frozen external velocity, solve the
inviscid problem for the transpiration it implies, relax, repeat. That loop is
what fails where the integral equations are singular, and both reference
implementations say so directly. Lokatt and Eller note that "the integral
boundary layer system is singular in separated flow regions" and that "the
singularity can be avoided by coupling the boundary layer equations to an
inviscid flow field model and solving both systems simultaneously", declaring
that out of their scope and restricting themselves to attached flow. Drela calls
retaining the simultaneous solution "highly desirable, if not essential" in 3D,
and reports that separation lines are then captured with no special treatment,
being where the characteristics of the hyperbolic system converge.

The joint unknown is the boundary-layer state of every panel stacked with the
source strengths,

```math
\mathbf{x}=\left(\mathbf{U},\;\boldsymbol\sigma\right),
\qquad
\mathbf{R}(\mathbf{x})=
\begin{pmatrix}
\mathbf{R}_v\!\left(\mathbf{U},\;\mathbf{q}_e(\boldsymbol\sigma)\right)\\[2pt]
\mathsf{K}\boldsymbol\sigma-\mathbf{v}_n(\mathbf{U})
\end{pmatrix},
```

with ``\mathbf{q}_e=\mathbf{q}_{rigid}+\mathsf{P}\mathsf{G}\boldsymbol\sigma``
the external velocity including the transpiration correction, tangentially
projected, and ``\mathbf{v}_n`` the transpiration velocity the layer implies.
The point is that the off-diagonal blocks are present: the layer sees the
external flow respond to its own displacement *within* the linear solve, which
is what relieves the singularity rather than iterating around it.

At boundary-element resolution the joint system is a few thousand unknowns, so
it is assembled and factorised densely. The inviscid influence is a constant
matrix (see [`indirect_gradient_operator`](@ref)), so no boundary-element solve
is repeated inside the iteration.
"""
function solve_simultaneous_coupling(grid::StructuredPanelGrid,
        rigid_edge_velocity::AbstractMatrix, kinematic_viscosity::Real;
        green_functions = (Rankine(), RankineReflected()), rho::Real = SETTINGS.rho,
        closure = ThreeDimensionalClosure(), reference = (0, 0, 0),
        newton_steps::Integer = 40, tolerance::Real = 1e-8,
        step_limit::Real = 1.0, source_scale = nothing,
        upwind_smoothing::Real = 0.15, kwargs...)
    mesh = grid.mesh
    size(rigid_edge_velocity) == (mesh.nfaces, 3) || throw(DimensionMismatch(
        "rigid_edge_velocity must have size (mesh.nfaces, 3)"))
    kinematic_viscosity > 0 ||
        throw(ArgumentError("kinematic_viscosity must be positive"))
    newton_steps >= 0 || throw(ArgumentError("newton_steps must be nonnegative"))

    wavenumber = zero(mesh.centers[1, 1])
    _, stiffness = assemble_matrices(green_functions, mesh, wavenumber; direct = false)
    influence = real.(stiffness)
    operator = indirect_gradient_operator(green_functions, mesh, wavenumber)

    topology = build_surface_topology(mesh)
    metrics = build_surface_metrics(mesh, topology)

    # Start from the uncoupled solve: the layer against the rigid velocity, and
    # the sources that its transpiration implies.
    seed = solve_surface_boundary_layer(mesh, rigid_edge_velocity,
        kinematic_viscosity; rho, reference, closure, topology, metrics,
        implicit_derivative = false, upwind_smoothing, kwargs...)
    sources = influence \ seed.transpiration_velocity
    states = copy(seed.diagnostics.states)
    scale = isnothing(source_scale) ?
            max(maximum(abs, sources), sqrt(eps(Float64))) : source_scale

    unknowns = vcat(vec(states), sources ./ scale)
    lagged = seed.shear_coefficient

    function joint_residual(vector)
        element = eltype(vector)
        panel_states = reshape(@view(vector[1:(3 * mesh.nfaces)]), 3, mesh.nfaces)
        strengths = @view(vector[(3 * mesh.nfaces + 1):end]) .* scale
        correction = _correction_velocity(operator, mesh, strengths)
        edge_velocity = rigid_edge_velocity .+ correction
        cache = build_surface_cache(mesh, edge_velocity, kinematic_viscosity; closure,
            topology, metrics, upwind_smoothing)
        viscous = global_residual(panel_states, mesh, cache,
            inflow_states(mesh, cache), lagged; couple_shear = true)
        transpiration = _transpiration_of(panel_states, mesh, cache, lagged)
        inviscid = influence * collect(element, strengths) .- transpiration
        return vcat(viscous, inviscid)
    end

    residual = joint_residual(unknowns)
    reference_norm = max(norm(residual), sqrt(eps(Float64)))
    completed = 0
    linear = Ref{Any}(nothing)
    for _ in 1:newton_steps
        norm(residual) / reference_norm < tolerance && break
        jacobian = ForwardDiff.jacobian(joint_residual, unknowns)
        step = _newton_step(linear, jacobian, residual)
        isnothing(step) && break
        accepted = false
        damping = one(eltype(step))
        for _ in 1:20
            # Clip the state part per panel, exactly as the uncoupled solve
            # does. A global admissibility test lets one panel veto the step for
            # every other, which on a hull collapses the damping to nothing and
            # stalls the iteration outright.
            candidate = copy(unknowns)
            for panel in 1:mesh.nfaces
                range = (3 * (panel - 1) + 1):(3 * panel)
                own = SVector{3}(@view unknowns[range])
                move = damping * SVector{3}(@view step[range])
                largest = maximum(abs, move)
                largest > step_limit && (move = move * (step_limit / largest))
                for _ in 1:20
                    _is_admissible(own + move, closure) && break
                    move = move / 2
                end
                candidate[range] .= own .+ move
            end
            sources_range = (3 * mesh.nfaces + 1):length(unknowns)
            candidate[sources_range] .= unknowns[sources_range] .+
                                        damping .* @view(step[sources_range])
            trial = joint_residual(candidate)
            if all(isfinite, trial) && norm(trial) < norm(residual)
                unknowns = candidate
                residual = trial
                accepted = true
                break
            end
            damping /= 2
        end
        completed += 1
        accepted || break
    end

    strengths = unknowns[(3 * mesh.nfaces + 1):end] .* scale
    edge_velocity = rigid_edge_velocity .+
                    _correction_velocity(operator, mesh, strengths)
    layer = solve_surface_boundary_layer(mesh, edge_velocity, kinematic_viscosity;
        rho, reference, closure, topology, metrics,
        initial = reshape(unknowns[1:(3 * mesh.nfaces)], 3, mesh.nfaces), kwargs...)
    relative = norm(residual) / reference_norm
    return SimultaneousCouplingResult(layer, strengths, edge_velocity, relative,
        completed, relative < tolerance)
end

# Transpiration velocity implied by a boundary-layer state, as the surface
# divergence of the mass defect with both components of the displacement
# thickness. Mirrors the assembly in `assemble_surface_result`, but takes the
# state directly so it can sit inside a differentiated joint residual.
function _transpiration_of(states, mesh::Mesh, cache::SurfaceBoundaryLayerCache,
        shear)
    element = promote_type(eltype(states), eltype(cache.speed))
    transpiration = zeros(element, mesh.nfaces)
    defects = Vector{SVector{2, element}}(undef, mesh.nfaces)
    for panel in 1:mesh.nfaces
        state = panel_fluxes(SVector{3}(@view states[:, panel]), cache.speed[panel],
            cache.cosine[panel], cache.sine[panel], shear[panel], cache.viscosity,
            cache.closure)
        defects[panel] = @SVector [cache.speed[panel] * state.displacement[1],
            cache.speed[panel] * state.displacement[2]]
    end
    for panel in 1:mesh.nfaces
        total = zero(element)
        for side in 1:4
            edge = cache.topology.cell_edges[panel, side]
            edge == 0 && continue
            donor = cache.edge_donor[edge]
            source = donor == 0 ? panel : donor
            defect = defects[source]
            if source != panel
                defect = rotation_into(cache.topology, cache.metrics, edge, panel) *
                         defect
            end
            area_vector = edge_area_vector_for(cache.topology, cache.metrics, edge,
                panel)
            total += defect[1] * area_vector[1] + defect[2] * area_vector[2]
        end
        transpiration[panel] = total / mesh.areas[panel]
    end
    return transpiration
end

solve_simultaneous_coupling(mesh::Mesh, args...; kwargs...) =
    throw(ArgumentError("solve_simultaneous_coupling needs a StructuredPanelGrid"))
