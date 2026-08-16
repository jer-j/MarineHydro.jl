raw"""
    VirtualCaptiveRun

Loads recorded at one combination of sway velocity and yaw rate.

`sway_force` and `yaw_moment` are the totals; `shear_force` and `shear_moment`
are the wall-stress part alone, so the difference is the part carried by the
displacement thickness through the transpiration potential. `separated` counts
the panels the closure reports as separated, and `converged` records whether the
layer solve met its own tolerance — a virtual test whose solves did not converge
tells you nothing, so it is kept alongside the loads rather than discarded.
"""
struct VirtualCaptiveRun{T}
    sway_velocity::T
    yaw_rate::T
    sway_force::T
    yaw_moment::T
    shear_force::T
    shear_moment::T
    separated::Int
    converged::Bool
end

raw"""
    VirtualCaptiveTest

A complete virtual captive-test matrix and the coefficients fitted to it.

`runs` holds every combination that was evaluated, `straight` the
zero-drift zero-rate reference that is subtracted from all of them, and
`coefficients` the regression described in [`fit_captive_coefficients`](@ref).
"""
struct VirtualCaptiveTest{T, C}
    runs::Vector{VirtualCaptiveRun{T}}
    straight::VirtualCaptiveRun{T}
    coefficients::C
    forward_speed::T
    length::T
end

raw"""
    virtual_captive_test(
        grid,
        surge_gradient,
        sway_gradient,
        yaw_gradient,
        forward_speed,
        kinematic_viscosity;
        drift_angles = range(-0.12, 0.12; length = 9),
        yaw_rates = nothing,
        rho = SETTINGS.rho,
        x_reference = 0,
        y_reference = 0,
        closure = HeadTurbulentClosure(),
        green_functions = (Rankine(), RankineReflected()),
        kwargs...
    )

Run a virtual captive test: sweep the hull through a matrix of drift angles and
yaw rates, record the sway force and yaw moment at each, and regress the
maneuvering coefficients from the result.

This is the numerical counterpart of an oblique-towing and rotating-arm
programme, and it answers something a linearisation cannot. `Y_v` and `N_v` as
computed by [`viscous_maneuvering_correction`](@ref) are derivatives *at* zero
drift; a captive test measures the load at drift angles a ship actually sees,
where the response is not linear. The fit therefore carries a cubic term
alongside the linear one, and the gap between the two is a direct statement of
how far the linear derivative can be trusted.

At each state the body-relative edge velocity is assembled exactly as the
correction assembles it, the boundary layer is solved, and its transpiration is
put through the same boundary-element solve, so that

```math
Y=\underbrace{\int_{S_h}\tau_y\,dS}_{\text{wall stress}}
-\rho U\int_{S_h}\frac{\partial\phi_Q}{\partial x}n_y\,dS,
```

with the pressure integral taken over the attached region only, following Wang
et al. The straight-ahead run is subtracted from every state, which is what a
towing tank does when it tares its balance.

The loads are the viscous contribution. Ideal flow at drift adds nothing to `Y`
and a Munk moment to `N`; combine with a potential-flow result through
[`apply_viscous_correction`](@ref) if the total is wanted.

`drift_angles` are in radians and give ``v=-U\sin\beta``. `yaw_rates` default to
a sweep of the same size in nondimensional ``r'=rL/U``; pass a single-element
range to hold one of them fixed.
"""
function virtual_captive_test(grid::StructuredPanelGrid, surge_gradient, sway_gradient,
        yaw_gradient, forward_speed::Real, kinematic_viscosity::Real;
        drift_angles = range(-0.12, 0.12; length = 9),
        yaw_rates = range(-0.3, 0.3; length = 9), rho::Real = SETTINGS.rho,
        x_reference::Real = 0, y_reference::Real = 0,
        closure = HeadTurbulentClosure(),
        green_functions = (Rankine(), RankineReflected()),
        minimum_edge_speed = nothing, integration_substeps::Integer = 8,
        stop_at_separation::Bool = true)
    mesh = grid.mesh
    forward_speed > 0 || throw(ArgumentError("forward_speed must be positive"))
    kinematic_viscosity > 0 ||
        throw(ArgumentError("kinematic_viscosity must be positive"))
    rho > 0 || throw(ArgumentError("rho must be positive"))
    isempty(drift_angles) && throw(ArgumentError("drift_angles must not be empty"))
    isempty(yaw_rates) && throw(ArgumentError("yaw_rates must not be empty"))

    ship_length = maximum(mesh.centers[:, 1]) - minimum(mesh.centers[:, 1])
    ship_length > 0 ||
        throw(ArgumentError("mesh must have nonzero longitudinal extent"))

    wavenumber = zero(mesh.centers[1, 1])
    S, K = assemble_matrices(green_functions, mesh, wavenumber; direct = false)
    geometry = _layer_geometry(closure, mesh)
    yaw_mode = maneuvering_boundary_condition(mesh, :yaw; x_reference, y_reference)

    warm_start = Ref{Any}(nothing)
    function evaluate(sway_velocity, yaw_rate)
        edge_velocity = body_relative_edge_velocity(mesh, forward_speed,
            sway_velocity, yaw_rate; surge_gradient, sway_gradient, yaw_gradient,
            x_reference, y_reference)
        layer = _solve_layer(closure, grid, edge_velocity, kinematic_viscosity; rho,
            reference = (x_reference, y_reference, zero(x_reference)),
            minimum_edge_speed, integration_substeps, stop_at_separation,
            geometry, initial = warm_start[])
        warm_start[] = _layer_warm_start(closure, layer)
        _, sources = solve(K, S, layer.transpiration_velocity; direct = false)
        gradient = real.(evaluate_indirect_potential_gradient(green_functions, mesh,
            sources, wavenumber))
        weight = mesh.areas .* .!layer.separated
        pressure_force = -rho * forward_speed *
                         dot(@view(gradient[:, 1]), mesh.normals[:, 2] .* weight)
        pressure_moment = -rho * forward_speed *
                          dot(@view(gradient[:, 1]), yaw_mode .* weight)
        return VirtualCaptiveRun(promote(sway_velocity, yaw_rate,
            layer.force[2] + pressure_force, layer.moment[3] + pressure_moment,
            layer.force[2], layer.moment[3])...,
            count(layer.separated), _layer_converged(layer))
    end

    zero_speed = zero(forward_speed)
    straight = evaluate(zero_speed, zero_speed)
    runs = [evaluate(-forward_speed * sin(angle), rate * forward_speed / ship_length)
            for angle in drift_angles, rate in yaw_rates]

    coefficients = fit_captive_coefficients(vec(runs), straight, forward_speed,
        ship_length; rho)
    return VirtualCaptiveTest(vec(runs), straight, coefficients,
        promote(forward_speed, ship_length)...)
end

# A march reports no diagnostics; treat it as converged, which it is by
# construction once it reaches the last station.
_layer_converged(result::ThreeDimensionalBoundaryLayerResult) =
    result.diagnostics.converged
_layer_converged(::Quasi3DBoundaryLayerResult) = true

raw"""
    fit_captive_coefficients(runs, straight, forward_speed, length; rho)

Regress maneuvering coefficients from a captive-test matrix.

Fits, in nondimensional form and after taring on the straight-ahead run,

```math
Y'=Y'_v v'+Y'_{vvv}v'^3+Y'_r r'+Y'_{rrr}r'^3+Y'_{vrr}v'r'^2+Y'_{vvr}v'^2r',
```

and the same for ``N'``, by ordinary least squares. The cubic and coupled terms
are the ones a linearisation about straight-ahead cannot reach, and they are
what the sweep exists to produce. Terms whose regressor is absent from the
matrix — a single-valued sweep in one variable, say — are dropped rather than
fitted to noise, so a pure drift sweep returns the sway family alone.

`residual` reports the root-mean-square misfit relative to the spread of the
data, so that a fit which the chosen basis cannot represent is visible instead
of silent.
"""
function fit_captive_coefficients(runs::AbstractVector{<:VirtualCaptiveRun},
        straight::VirtualCaptiveRun, forward_speed::Real, length::Real;
        rho::Real = SETTINGS.rho)
    sway = [run.sway_velocity / forward_speed for run in runs]
    yaw = [run.yaw_rate * length / forward_speed for run in runs]
    force_scale = rho / 2 * length^2 * forward_speed^2
    moment_scale = rho / 2 * length^3 * forward_speed^2
    force = [(run.sway_force - straight.sway_force) / force_scale for run in runs]
    moment = [(run.yaw_moment - straight.yaw_moment) / moment_scale for run in runs]

    candidates = (:v => sway, :vvv => sway .^ 3, :r => yaw, :rrr => yaw .^ 3,
        :vrr => sway .* yaw .^ 2, :vvr => sway .^ 2 .* yaw)
    active = [(name, column) for (name, column) in candidates
              if maximum(abs, column) > sqrt(eps(Float64)) * maximum(abs, sway .+ yaw)]
    isempty(active) && throw(ArgumentError("captive matrix varies in neither v nor r"))
    design = reduce(hcat, [column for (_, column) in active])

    sway_fit = design \ force
    yaw_fit = design \ moment
    names = [name for (name, _) in active]
    spread(values) = maximum(values) - minimum(values)
    return (sway = NamedTuple{Tuple(Symbol.("Y_", names))}(Tuple(sway_fit)),
        yaw = NamedTuple{Tuple(Symbol.("N_", names))}(Tuple(yaw_fit)),
        residual = (sway = norm(design * sway_fit .- force) /
                           max(sqrt(Base.length(force)) * spread(force), eps()),
            yaw = norm(design * yaw_fit .- moment) /
                  max(sqrt(Base.length(moment)) * spread(moment), eps())))
end
