raw"""
    WangHydrodynamicDerivatives

Dimensional linear sway and yaw hydrodynamic derivatives using the sign and
normalization conventions of Wang, Bakountouzis, and Katory (2000).

The fields correspond to

```math
Y = Y_{\dot v}\dot v + Y_v v + Y_{\dot r}\dot r + Y_r r,
\qquad
N = N_{\dot v}\dot v + N_v v + N_{\dot r}\dot r + N_r r.
```

The signs follow Wang et al.: the panel normals point from the body into the
fluid, ``x`` points forward, and ``y`` points to port.
"""
struct WangHydrodynamicDerivatives{T}
    Y_vdot::T
    Y_v::T
    Y_rdot::T
    Y_r::T
    N_vdot::T
    N_v::T
    N_rdot::T
    N_r::T
end

raw"""
    MMGLinearHullDerivatives

Linear dimensionless sway-yaw hull derivatives in the standard MMG
normalization,

```math
Y_H'=Y_v'v'+Y_r'r',\qquad
N_H'=N_v'v'+N_r'r',
```

where ``v'=v/U`` and ``r'=rL/U``. The standard force and moment scales are
``\tfrac12\rho LTU^2`` and ``\tfrac12\rho L^2TU^2``, respectively.
"""
struct MMGLinearHullDerivatives{T}
    Y_v::T
    Y_r::T
    N_v::T
    N_r::T
end

function MMGLinearHullDerivatives(Y_v, Y_r, N_v, N_r)
    values = promote(Y_v, Y_r, N_v, N_r)
    return MMGLinearHullDerivatives{typeof(first(values))}(values...)
end

raw"""
    mmg_to_wang_velocity_derivatives(derivatives, draft, length)

Convert standard-MMG linear hull derivatives to Wang's dimensionless scales.
With ``\lambda_T=T/L``, all four first-order velocity derivatives satisfy

```math
\left(Y_v',Y_r',N_v',N_r'\right)_W
=\lambda_T\left(Y_v',Y_r',N_v',N_r'\right)_{MMG}.
```

This changes only the normalization. Coordinate-sign conversions, if needed,
must be handled before calling this function.
"""
function mmg_to_wang_velocity_derivatives(
    derivatives::MMGLinearHullDerivatives,
    draft::Real,
    length::Real,
)
    draft > 0 || throw(ArgumentError("draft must be positive"))
    length > 0 || throw(ArgumentError("length must be positive"))
    draft_ratio = draft / length
    return MMGLinearHullDerivatives(
        draft_ratio * derivatives.Y_v,
        draft_ratio * derivatives.Y_r,
        draft_ratio * derivatives.N_v,
        draft_ratio * derivatives.N_r,
    )
end

function WangHydrodynamicDerivatives(
    Y_vdot,
    Y_v,
    Y_rdot,
    Y_r,
    N_vdot,
    N_v,
    N_rdot,
    N_r,
)
    values = promote(Y_vdot, Y_v, Y_rdot, Y_r, N_vdot, N_v, N_rdot, N_r)
    return WangHydrodynamicDerivatives{typeof(first(values))}(values...)
end

raw"""
    WangManeuveringResult

Solution data and dimensional derivatives from [`solve_wang_maneuvering`](@ref).

`boundary_residual` is the relative residual of the indirect boundary-integral
system. `velocity_mask` identifies the panels used for the velocity-dependent
derivatives. Acceleration-dependent derivatives always use the complete mesh.
"""
struct WangManeuveringResult{P,S,G,D,M,R}
    potential::P
    sources::S
    potential_gradient::G
    derivatives::D
    velocity_mask::M
    boundary_residual::R
end

raw"""
    PotentialFlowManeuveringResult

Solution data from [`solve_potential_flow_maneuvering`](@ref), including the
independent unit-sway and unit-yaw potentials.

`boundary_residuals.sway` and `boundary_residuals.yaw` are the relative
residuals of the two indirect boundary-integral systems. `velocity_mask`
identifies panels used for velocity-dependent derivatives. All
acceleration-dependent derivatives use the complete mesh.
"""
struct PotentialFlowManeuveringResult{P,S,G,D,M,R}
    potentials::P
    sources::S
    potential_gradients::G
    derivatives::D
    velocity_mask::M
    boundary_residuals::R
end

raw"""
    RigidBodyPotentialResult

The unit-potential solution for one rigid-body mode. `potential_gradient` is
the fluid-side surface trace, and `boundary_residual` is the relative residual
of the indirect boundary-integral system.
"""
struct RigidBodyPotentialResult{P,S,G,R}
    potential::P
    sources::S
    potential_gradient::G
    boundary_residual::R
end

raw"""
    WangRestrictedWaterResult

Sectionwise restricted-water correction from Wang et al. (2000).
`local_speed` is the accelerated longitudinal speed ``U_1(x)``, `elevation`
is the local mean free-surface elevation ``\zeta(x)``, and `mean_elevation`
is its waterline-beam-weighted longitudinal mean.
"""
struct WangRestrictedWaterResult{V,T}
    local_speed::V
    elevation::V
    mean_elevation::T
end

function _trapezoidal_integral(coordinate, values)
    integral_value = zero(promote_type(eltype(coordinate), eltype(values)))
    for station in 1:(length(coordinate) - 1)
        interval = abs(coordinate[station + 1] - coordinate[station])
        integral_value += interval *
            (values[station] + values[station + 1]) / 2
    end
    return integral_value
end

raw"""
    wang_restricted_water_elevation(
        coordinate,
        section_area,
        waterline_beam,
        forward_speed;
        channel_area=Inf,
        channel_surface_width=Inf,
        gravity=SETTINGS.g,
    )

Solve Wang et al.'s sectionwise continuity and Bernoulli correction,

```math
U_0A_0=U_1(x)\left[A_0-S(x)+\zeta(x)
\left(W-B(x)\right)\right],
\qquad
\frac{U_0^2}{2}=\frac{U_1(x)^2}{2}+g\zeta(x).
```

Here ``A_0`` is the undisturbed channel cross-sectional area and ``W`` is its
surface width. The subcritical root nearest ``U_0`` is selected. Setting both
channel dimensions to `Inf` returns the unrestricted-water limit
``U_1=U_0`` and ``\zeta=0``.

This is the mean blockage and squat correction in section 2.2 of Wang et al.,
not a finite-Froude-number wave-making solution.
"""
function wang_restricted_water_elevation(
    coordinate::AbstractVector,
    section_area::AbstractVector,
    waterline_beam::AbstractVector,
    forward_speed::Real;
    channel_area::Real=Inf,
    channel_surface_width::Real=Inf,
    gravity::Real=SETTINGS.g,
)
    station_count = length(coordinate)
    station_count >= 2 || throw(ArgumentError(
        "at least two longitudinal stations are required",
    ))
    length(section_area) == station_count || throw(DimensionMismatch(
        "coordinate and section_area must have equal length",
    ))
    length(waterline_beam) == station_count || throw(DimensionMismatch(
        "coordinate and waterline_beam must have equal length",
    ))
    coordinate_difference = diff(coordinate)
    all(coordinate_difference .> 0) || all(coordinate_difference .< 0) ||
        throw(ArgumentError("coordinate must be strictly monotone"))
    all(section_area .>= 0) || throw(ArgumentError(
        "section_area must be nonnegative",
    ))
    all(waterline_beam .>= 0) || throw(ArgumentError(
        "waterline_beam must be nonnegative",
    ))
    forward_speed > 0 || throw(ArgumentError(
        "forward_speed must be positive",
    ))
    gravity > 0 || throw(ArgumentError("gravity must be positive"))

    unrestricted = isinf(channel_area) && isinf(channel_surface_width)
    isinf(channel_area) == isinf(channel_surface_width) || throw(ArgumentError(
        "channel_area and channel_surface_width must both be finite or both Inf",
    ))
    T = promote_type(
        eltype(coordinate),
        eltype(section_area),
        eltype(waterline_beam),
        typeof(forward_speed),
        typeof(gravity),
    )
    if unrestricted
        local_speed = fill(convert(T, forward_speed), station_count)
        elevation = zeros(T, station_count)
        return WangRestrictedWaterResult(local_speed, elevation, zero(T))
    end

    channel_area > maximum(section_area) || throw(ArgumentError(
        "channel_area must exceed every immersed section area",
    ))
    channel_surface_width > maximum(waterline_beam) || throw(ArgumentError(
        "channel_surface_width must exceed every waterline beam",
    ))
    local_speed = Vector{T}(undef, station_count)
    elevation = Vector{T}(undef, station_count)
    for station in eachindex(coordinate)
        area = section_area[station]
        beam = waterline_beam[station]
        if iszero(area)
            local_speed[station] = forward_speed
            elevation[station] = zero(T)
            continue
        end
        width_difference = channel_surface_width - beam
        surface_elevation(speed) =
            (forward_speed^2 - speed^2) / (2 * gravity)
        residual(speed) = speed * (
            channel_area - area +
            surface_elevation(speed) * width_difference
        ) - forward_speed * channel_area
        residual_derivative(speed) =
            channel_area - area + surface_elevation(speed) * width_difference -
            speed^2 * width_difference / gravity

        lower_speed = convert(T, forward_speed)
        upper_speed = lower_speed * (one(T) + convert(T, 1e-3))
        bracketed = false
        for _ in 1:100
            if residual(upper_speed) >= 0
                bracketed = true
                break
            end
            residual_derivative(upper_speed) > 0 || break
            upper_speed = forward_speed + 1.5 *
                (upper_speed - forward_speed)
        end
        bracketed || throw(DomainError(
            area,
            "no subcritical restricted-water root exists at this section",
        ))
        for _ in 1:80
            midpoint = (lower_speed + upper_speed) / 2
            if residual(midpoint) < 0
                lower_speed = midpoint
            else
                upper_speed = midpoint
            end
        end
        local_speed[station] = (lower_speed + upper_speed) / 2
        elevation[station] = surface_elevation(local_speed[station])
    end
    beam_integral = _trapezoidal_integral(coordinate, waterline_beam)
    mean_elevation = iszero(beam_integral) ? zero(T) :
        _trapezoidal_integral(
            coordinate,
            elevation .* waterline_beam,
        ) / beam_integral
    return WangRestrictedWaterResult(
        local_speed,
        elevation,
        mean_elevation,
    )
end

raw"""
    maneuvering_boundary_condition(
        mesh,
        motion;
        x_reference=0,
        y_reference=0,
    )

Return a unit rigid-body Neumann boundary condition for surge, sway, or yaw.
For `motion=:surge` and `motion=:sway`, respectively,

```math
g_U = n_x,
\qquad
g_v = n_y.
```

For `motion=:yaw`, the exact three-dimensional mode about
``(x_0,y_0)`` is

```math
g_r = -(y-y_0)n_x + (x-x_0)n_y.
```

The mesh normals must point from the body into the fluid. Coordinates follow
the Wang convention: ``x`` points forward, ``y`` points to port, and positive
yaw is about the positive ``z`` axis.
"""
function maneuvering_boundary_condition(
    mesh::Mesh,
    motion::Symbol;
    x_reference::Real=0,
    y_reference::Real=0,
)
    if motion === :surge
        return vec(mesh.normals[:, 1])
    elseif motion === :sway
        return vec(mesh.normals[:, 2])
    elseif motion === :yaw
        x = mesh.centers[:, 1] .- x_reference
        y = mesh.centers[:, 2] .- y_reference
        nx = mesh.normals[:, 1]
        ny = mesh.normals[:, 2]
        return -y .* nx .+ x .* ny
    end
    throw(ArgumentError("motion must be :surge, :sway, or :yaw"))
end

maneuvering_boundary_condition(floatingbody::FloatingBody, args...; kwargs...) =
    maneuvering_boundary_condition(floatingbody.mesh, args...; kwargs...)

raw"""
    wang_sway_boundary_condition(mesh)

Return the unit-sway Neumann boundary condition

```math
\frac{\partial \phi_2}{\partial n} = n_y.
```

The mesh normals must point from the body into the fluid.
"""
wang_sway_boundary_condition(mesh::Mesh) =
    maneuvering_boundary_condition(mesh, :sway)

raw"""
    wang_stern_mask(mesh, x_cut; x_reference=0, bow_direction=:positive_x)

Construct the panel mask used for the Schmitz stern truncation in the Wang
velocity-derivative integrals.

`x_cut` is measured relative to `x_reference`. With
`bow_direction=:positive_x`, panels satisfying ``x-x_0 \ge x_{cut}`` are kept.
Use `bow_direction=:negative_x` when the mesh longitudinal axis points aft.
"""
function wang_stern_mask(
    mesh::Mesh,
    x_cut::Real;
    x_reference::Real=0,
    bow_direction::Symbol=:positive_x,
)
    x = mesh.centers[:, 1] .- x_reference
    if bow_direction === :positive_x
        return BitVector(x .>= x_cut)
    elseif bow_direction === :negative_x
        return BitVector(x .<= x_cut)
    end
    throw(ArgumentError("bow_direction must be :positive_x or :negative_x"))
end

wang_stern_mask(floatingbody::FloatingBody, args...; kwargs...) =
    wang_stern_mask(floatingbody.mesh, args...; kwargs...)

raw"""
    evaluate_indirect_potential_gradient(
        green_functions,
        mesh,
        sources,
        wavenumber=0,
    )

Evaluate the fluid-side trace of the potential gradient generated by an
indirect single-layer source distribution.

For a target panel ``i``, the calculation is

```math
\nabla\phi_i = c_i\boldsymbol{n}_i\sigma_i
- \frac{1}{4\pi}\sum_j
\int_{S_j}\nabla_{\boldsymbol{x}}G(\boldsymbol{x}_i,\boldsymbol{\xi})
\sigma_j\,\mathrm{d}S_{\xi},
```

using the same jump convention as [`assemble_matrices`](@ref). The result is an
``N \times 3`` matrix whose columns are ``\phi_x``, ``\phi_y``, and ``\phi_z``.
"""
function evaluate_indirect_potential_gradient(
    green_functions,
    mesh::Mesh,
    sources::AbstractVector,
    wavenumber::Real=0,
)
    length(sources) == mesh.nfaces || throw(DimensionMismatch(
        "sources must contain one value per mesh face",
    ))

    free_surface = zero(mesh.centers[1, 3])
    gradient_vectors = [
        let
            target = element(mesh, i)
            one_coordinate = one(mesh.centers[i, 3])
            jump = abs(mesh.centers[i, 3] - free_surface) < 1e-8 ?
                one_coordinate : one_coordinate / 2
            principal_value = sum(
                integral_gradient(
                    green_functions,
                    target,
                    element(mesh, j),
                    wavenumber;
                    with_respect_to_first_variable=true,
                ) .* sources[j]
                for j in 1:mesh.nfaces
            )
            jump .* normal(target) .* sources[i] .- principal_value ./ (2τ̅)
        end
        for i in 1:mesh.nfaces
    ]

    return copy(stack(gradient_vectors)')
end

raw"""
    solve_rigid_body_potential(
        mesh,
        motion;
        x_reference=0,
        y_reference=0,
        green_functions=(Rankine(), RankineReflected()),
    )

Solve a unit surge, sway, or yaw potential on the unrestricted double-body
domain. This lower-level entry point is useful when an external-flow field is
needed in addition to integrated maneuvering derivatives, for example as the
edge condition for a boundary-layer calculation.
"""
function solve_rigid_body_potential(
    mesh::Mesh,
    motion::Symbol;
    x_reference::Real=0,
    y_reference::Real=0,
    green_functions=(Rankine(), RankineReflected()),
)
    mesh.nfaces > 0 || throw(ArgumentError("mesh must contain at least one face"))
    boundary_condition = maneuvering_boundary_condition(
        mesh,
        motion;
        x_reference,
        y_reference,
    )
    wavenumber = zero(mesh.centers[1, 1])
    S, K = assemble_matrices(
        green_functions,
        mesh,
        wavenumber;
        direct=false,
    )
    potential, sources = solve(K, S, boundary_condition; direct=false)
    potential_gradient = evaluate_indirect_potential_gradient(
        green_functions,
        mesh,
        sources,
        wavenumber,
    )

    residual_denominator = norm(boundary_condition)
    residual_norm = norm(K * sources - boundary_condition)
    boundary_residual = iszero(residual_denominator) ? residual_norm :
        residual_norm / residual_denominator
    return RigidBodyPotentialResult(
        real.(potential),
        real.(sources),
        real.(potential_gradient),
        boundary_residual,
    )
end

solve_rigid_body_potential(floatingbody::FloatingBody, args...; kwargs...) =
    solve_rigid_body_potential(floatingbody.mesh, args...; kwargs...)

raw"""
    wang_hydrodynamic_derivatives(
        mesh,
        potential,
        longitudinal_gradient,
        forward_speed;
        rho=SETTINGS.rho,
        x_reference=0,
        velocity_mask=nothing,
    )

Evaluate the dimensional Wang sway-yaw derivatives from the unit-sway
potential ``\phi_2`` and its longitudinal derivative ``\partial\phi_2/\partial
x``.

`velocity_mask` implements the semi-empirical Schmitz stern truncation. Pass
`nothing` to integrate the velocity-dependent derivatives over the whole hull.
The acceleration-dependent derivatives always use the whole hull.
"""
function wang_hydrodynamic_derivatives(
    mesh::Mesh,
    potential::AbstractVector,
    longitudinal_gradient::AbstractVector,
    forward_speed::Real;
    rho::Real=SETTINGS.rho,
    x_reference::Real=0,
    velocity_mask::Union{Nothing,AbstractVector{Bool}}=nothing,
)
    length(potential) == mesh.nfaces || throw(DimensionMismatch(
        "potential must contain one value per mesh face",
    ))
    length(longitudinal_gradient) == mesh.nfaces || throw(DimensionMismatch(
        "longitudinal_gradient must contain one value per mesh face",
    ))
    forward_speed >= 0 || throw(ArgumentError("forward_speed must be nonnegative"))
    rho > 0 || throw(ArgumentError("rho must be positive"))

    mask = isnothing(velocity_mask) ? trues(mesh.nfaces) : BitVector(velocity_mask)
    length(mask) == mesh.nfaces || throw(DimensionMismatch(
        "velocity_mask must contain one value per mesh face",
    ))

    x = mesh.centers[:, 1] .- x_reference
    ny = mesh.normals[:, 2]
    dS = mesh.areas
    phi = potential
    phi_x = longitudinal_gradient
    xphi_x = phi .+ x .* phi_x

    acceleration_weight = ny .* dS
    velocity_weight = acceleration_weight .* mask

    Y_vdot = rho * sum(phi .* acceleration_weight)
    Y_rdot = rho * sum(x .* phi .* acceleration_weight)
    N_vdot = Y_rdot
    N_rdot = rho * sum(x .^ 2 .* phi .* acceleration_weight)

    Y_v = -rho * forward_speed * sum(phi_x .* velocity_weight)
    Y_r = -rho * forward_speed * sum(xphi_x .* velocity_weight)
    N_v = -rho * forward_speed * sum(x .* phi_x .* velocity_weight)
    N_r = -rho * forward_speed * sum(x .* xphi_x .* velocity_weight)

    return WangHydrodynamicDerivatives(
        Y_vdot,
        Y_v,
        Y_rdot,
        Y_r,
        N_vdot,
        N_v,
        N_rdot,
        N_r,
    )
end

raw"""
    potential_flow_hydrodynamic_derivatives(
        mesh,
        sway_potential,
        yaw_potential,
        sway_longitudinal_gradient,
        yaw_longitudinal_gradient,
        forward_speed;
        rho=SETTINGS.rho,
        x_reference=0,
        y_reference=0,
        velocity_mask=nothing,
    )

Evaluate the dimensional linear sway-yaw derivatives from independently
solved unit-sway and unit-yaw potentials. With

```math
g_v=n_y,
\qquad
g_r=-(y-y_0)n_x+(x-x_0)n_y,
```

the acceleration-dependent coefficients are

```math
Q_{\dot\jmath}=\rho\int_{S_h}\phi_j g_Q\,\mathrm{d}S,
```

and the zero-frequency convective coefficients are

```math
Q_j=-\rho U\int_{S_h^*}\phi_{j,x}g_Q\,\mathrm{d}S.
```

Here ``j`` is sway or yaw and ``Q`` is sway force or yaw moment. This is the
exact rigid-body yaw extension of the strict Wang approximation
``\phi_r\approx(x-x_0)\phi_v``. `velocity_mask` affects only the convective
integrals.
"""
function potential_flow_hydrodynamic_derivatives(
    mesh::Mesh,
    sway_potential::AbstractVector,
    yaw_potential::AbstractVector,
    sway_longitudinal_gradient::AbstractVector,
    yaw_longitudinal_gradient::AbstractVector,
    forward_speed::Real;
    rho::Real=SETTINGS.rho,
    x_reference::Real=0,
    y_reference::Real=0,
    velocity_mask::Union{Nothing,AbstractVector{Bool}}=nothing,
)
    fields = (
        sway_potential=sway_potential,
        yaw_potential=yaw_potential,
        sway_longitudinal_gradient=sway_longitudinal_gradient,
        yaw_longitudinal_gradient=yaw_longitudinal_gradient,
    )
    for (name, values) in pairs(fields)
        length(values) == mesh.nfaces || throw(DimensionMismatch(
            "$(name) must contain one value per mesh face",
        ))
    end
    forward_speed >= 0 || throw(ArgumentError("forward_speed must be nonnegative"))
    rho > 0 || throw(ArgumentError("rho must be positive"))

    mask = isnothing(velocity_mask) ? trues(mesh.nfaces) : BitVector(velocity_mask)
    length(mask) == mesh.nfaces || throw(DimensionMismatch(
        "velocity_mask must contain one value per mesh face",
    ))

    g_v = maneuvering_boundary_condition(mesh, :sway)
    g_r = maneuvering_boundary_condition(
        mesh,
        :yaw;
        x_reference,
        y_reference,
    )
    dS = mesh.areas
    acceleration_sway_weight = g_v .* dS
    acceleration_yaw_weight = g_r .* dS
    velocity_sway_weight = acceleration_sway_weight .* mask
    velocity_yaw_weight = acceleration_yaw_weight .* mask

    Y_vdot = rho * sum(sway_potential .* acceleration_sway_weight)
    Y_rdot = rho * sum(yaw_potential .* acceleration_sway_weight)
    N_vdot = rho * sum(sway_potential .* acceleration_yaw_weight)
    N_rdot = rho * sum(yaw_potential .* acceleration_yaw_weight)

    Y_v = -rho * forward_speed * sum(
        sway_longitudinal_gradient .* velocity_sway_weight,
    )
    Y_r = -rho * forward_speed * sum(
        yaw_longitudinal_gradient .* velocity_sway_weight,
    )
    N_v = -rho * forward_speed * sum(
        sway_longitudinal_gradient .* velocity_yaw_weight,
    )
    N_r = -rho * forward_speed * sum(
        yaw_longitudinal_gradient .* velocity_yaw_weight,
    )

    return WangHydrodynamicDerivatives(
        Y_vdot,
        Y_v,
        Y_rdot,
        Y_r,
        N_vdot,
        N_v,
        N_rdot,
        N_r,
    )
end

raw"""
    nondimensionalize_wang_derivatives(derivatives, length, forward_speed; rho=SETTINGS.rho)

Apply the nondimensionalization used by Wang et al. (2000):

```math
Y'_{\dot v}=\frac{Y_{\dot v}}{\tfrac12\rho L^3},\quad
Y'_v=\frac{Y_v}{\tfrac12\rho L^2U},\quad
Y'_{\dot r}=\frac{Y_{\dot r}}{\tfrac12\rho L^4},\quad
Y'_r=\frac{Y_r}{\tfrac12\rho L^3U},
```

with analogous moment scalings for ``N``.
"""
function nondimensionalize_wang_derivatives(
    derivatives::WangHydrodynamicDerivatives,
    length::Real,
    forward_speed::Real;
    rho::Real=SETTINGS.rho,
)
    length > 0 || throw(ArgumentError("length must be positive"))
    forward_speed > 0 || throw(ArgumentError(
        "forward_speed must be positive for velocity-derivative normalization",
    ))
    rho > 0 || throw(ArgumentError("rho must be positive"))

    half_rho = rho / 2
    return WangHydrodynamicDerivatives(
        derivatives.Y_vdot / (half_rho * length^3),
        derivatives.Y_v / (half_rho * length^2 * forward_speed),
        derivatives.Y_rdot / (half_rho * length^4),
        derivatives.Y_r / (half_rho * length^3 * forward_speed),
        derivatives.N_vdot / (half_rho * length^4),
        derivatives.N_v / (half_rho * length^3 * forward_speed),
        derivatives.N_rdot / (half_rho * length^5),
        derivatives.N_r / (half_rho * length^4 * forward_speed),
    )
end

raw"""
    nondimensionalize_maneuvering_derivatives(args...; kwargs...)

Apply Wang's nondimensionalization to either strict-Wang or independent-yaw
derivatives. This is the method-neutral alias of
[`nondimensionalize_wang_derivatives`](@ref).
"""
nondimensionalize_maneuvering_derivatives(args...; kwargs...) =
    nondimensionalize_wang_derivatives(args...; kwargs...)

raw"""
    solve_wang_maneuvering(
        mesh,
        forward_speed;
        rho=SETTINGS.rho,
        x_reference=0,
        velocity_mask=nothing,
        green_functions=(Rankine(), RankineReflected()),
    )

Solve the unrestricted, infinite-depth, zero-frequency Wang maneuvering
problem and return the dimensional sway-yaw derivatives.

The default Green function is the rigid-free-surface, reflected-Rankine kernel

```math
G_0 = \frac{1}{r}+\frac{1}{r'},
```

in MarineHydro's internal normalization. Forward speed enters only through
the convective Bernoulli term used in the velocity derivatives; it is not used
in the boundary-integral operator. This is Wang's double-body approximation,
not a steady wave-making calculation.
"""
function solve_wang_maneuvering(
    mesh::Mesh,
    forward_speed::Real;
    rho::Real=SETTINGS.rho,
    x_reference::Real=0,
    velocity_mask::Union{Nothing,AbstractVector{Bool}}=nothing,
    green_functions=(Rankine(), RankineReflected()),
)
    forward_speed >= 0 || throw(ArgumentError("forward_speed must be nonnegative"))
    rho > 0 || throw(ArgumentError("rho must be positive"))
    mesh.nfaces > 0 || throw(ArgumentError("mesh must contain at least one face"))

    mask = isnothing(velocity_mask) ? trues(mesh.nfaces) : BitVector(velocity_mask)
    length(mask) == mesh.nfaces || throw(DimensionMismatch(
        "velocity_mask must contain one value per mesh face",
    ))

    boundary_condition = wang_sway_boundary_condition(mesh)
    # Forward speed is absent from Wang's zero-frequency boundary-value
    # problem. Keeping it out of the operator also avoids unnecessary dual
    # propagation when differentiating the velocity derivatives with respect
    # to speed.
    wavenumber = zero(mesh.centers[1, 1])
    S, K = assemble_matrices(
        green_functions,
        mesh,
        wavenumber;
        direct=false,
    )
    potential, sources = solve(K, S, boundary_condition; direct=false)
    potential_gradient = evaluate_indirect_potential_gradient(
        green_functions,
        mesh,
        sources,
        wavenumber,
    )

    real_potential = real.(potential)
    real_sources = real.(sources)
    real_gradient = real.(potential_gradient)
    derivatives = wang_hydrodynamic_derivatives(
        mesh,
        real_potential,
        real_gradient[:, 1],
        forward_speed;
        rho,
        x_reference,
        velocity_mask=mask,
    )

    residual_denominator = norm(boundary_condition)
    residual_norm = norm(K * sources - boundary_condition)
    boundary_residual = iszero(residual_denominator) ? residual_norm :
        residual_norm / residual_denominator

    return WangManeuveringResult(
        real_potential,
        real_sources,
        real_gradient,
        derivatives,
        mask,
        boundary_residual,
    )
end

raw"""
    solve_potential_flow_maneuvering(
        mesh,
        forward_speed;
        rho=SETTINGS.rho,
        x_reference=0,
        y_reference=0,
        velocity_mask=nothing,
        green_functions=(Rankine(), RankineReflected()),
    )

Solve independent unit-sway and unit-yaw zero-frequency boundary-value
problems on an unrestricted, infinite-depth domain. The two right-hand sides
share one assembled boundary-element operator.

The default reflected-Rankine kernel imposes a rigid free surface. Forward
speed enters only through the convective Bernoulli terms and is absent from
the boundary-integral operator. This method therefore remains a double-body
approximation rather than a finite-Froude-number wave-making calculation.
"""
function solve_potential_flow_maneuvering(
    mesh::Mesh,
    forward_speed::Real;
    rho::Real=SETTINGS.rho,
    x_reference::Real=0,
    y_reference::Real=0,
    velocity_mask::Union{Nothing,AbstractVector{Bool}}=nothing,
    green_functions=(Rankine(), RankineReflected()),
)
    forward_speed >= 0 || throw(ArgumentError("forward_speed must be nonnegative"))
    rho > 0 || throw(ArgumentError("rho must be positive"))
    mesh.nfaces > 0 || throw(ArgumentError("mesh must contain at least one face"))

    mask = isnothing(velocity_mask) ? trues(mesh.nfaces) : BitVector(velocity_mask)
    length(mask) == mesh.nfaces || throw(DimensionMismatch(
        "velocity_mask must contain one value per mesh face",
    ))

    sway_boundary_condition = maneuvering_boundary_condition(mesh, :sway)
    yaw_boundary_condition = maneuvering_boundary_condition(
        mesh,
        :yaw;
        x_reference,
        y_reference,
    )
    boundary_conditions = hcat(sway_boundary_condition, yaw_boundary_condition)
    wavenumber = zero(mesh.centers[1, 1])
    S, K = assemble_matrices(
        green_functions,
        mesh,
        wavenumber;
        direct=false,
    )
    potentials, sources = solve(K, S, boundary_conditions; direct=false)
    sway_gradient = evaluate_indirect_potential_gradient(
        green_functions,
        mesh,
        sources[:, 1],
        wavenumber,
    )
    yaw_gradient = evaluate_indirect_potential_gradient(
        green_functions,
        mesh,
        sources[:, 2],
        wavenumber,
    )

    real_potentials = real.(potentials)
    real_sources = real.(sources)
    real_sway_gradient = real.(sway_gradient)
    real_yaw_gradient = real.(yaw_gradient)
    derivatives = potential_flow_hydrodynamic_derivatives(
        mesh,
        real_potentials[:, 1],
        real_potentials[:, 2],
        real_sway_gradient[:, 1],
        real_yaw_gradient[:, 1],
        forward_speed;
        rho,
        x_reference,
        y_reference,
        velocity_mask=mask,
    )

    function relative_residual(column)
        denominator = norm(boundary_conditions[:, column])
        residual = norm(K * sources[:, column] - boundary_conditions[:, column])
        return iszero(denominator) ? residual : residual / denominator
    end

    return PotentialFlowManeuveringResult(
        (sway=real_potentials[:, 1], yaw=real_potentials[:, 2]),
        (sway=real_sources[:, 1], yaw=real_sources[:, 2]),
        (sway=real_sway_gradient, yaw=real_yaw_gradient),
        derivatives,
        mask,
        (sway=relative_residual(1), yaw=relative_residual(2)),
    )
end


solve_wang_maneuvering(floatingbody::FloatingBody, args...; kwargs...) =
    solve_wang_maneuvering(floatingbody.mesh, args...; kwargs...)

solve_potential_flow_maneuvering(floatingbody::FloatingBody, args...; kwargs...) =
    solve_potential_flow_maneuvering(floatingbody.mesh, args...; kwargs...)
