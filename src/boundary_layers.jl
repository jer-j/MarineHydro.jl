using LinearAlgebra: cross, dot, norm

raw"""
    HeadTurbulentClosure(; initial_shape_factor=1.4,
                          separation_shape_factor=2.4)

Parameters for Head's entrainment formulation of an attached turbulent
integral boundary layer. The state consists of momentum thickness ``\theta``
and displacement shape factor ``H=\delta^*/\theta``. A station is flagged as
separated when `shape_factor >= separation_shape_factor`. The default solver
freezes the integral state and sets wall shear to zero after that point because
the attached-flow closure is no longer valid.
"""
struct HeadTurbulentClosure{T}
    initial_shape_factor::T
    separation_shape_factor::T
end

function HeadTurbulentClosure(;
    initial_shape_factor::Real=1.4,
    separation_shape_factor::Real=2.4,
)
    initial_shape_factor > 1.1 || throw(ArgumentError(
        "initial_shape_factor must exceed 1.1 for the Head closure",
    ))
    separation_shape_factor > initial_shape_factor || throw(ArgumentError(
        "separation_shape_factor must exceed initial_shape_factor",
    ))
    values = promote(initial_shape_factor, separation_shape_factor)
    return HeadTurbulentClosure{typeof(first(values))}(values...)
end

raw"""
    IntegralBoundaryLayerResult

Stationwise solution of Head's turbulent integral boundary-layer equations.
`kinetic_shape_factor` stores Head's entrainment variable ``H_1`` and
`separated` identifies stations beyond the supplied attached-flow threshold.
"""
struct IntegralBoundaryLayerResult{S,U,T,H,K,D,C,B}
    coordinate::S
    edge_speed::U
    momentum_thickness::T
    shape_factor::H
    kinetic_shape_factor::K
    displacement_thickness::D
    skin_friction_coefficient::C
    separated::B
end

raw"""
    Quasi3DBoundaryLayerResult

One-way, strip-marched integral boundary-layer solution on a three-dimensional
panel surface. Scalar integral states are stored per panel. `wall_shear` is
the traction exerted by the fluid on the body and
`transpiration_velocity` is positive along the outward body normal.
"""
struct Quasi3DBoundaryLayerResult{V,M,B,F}
    edge_speed::V
    momentum_thickness::V
    displacement_thickness::V
    shape_factor::V
    skin_friction_coefficient::V
    wall_shear::M
    transpiration_velocity::V
    separated::B
    force::F
    moment::F
end

raw"""
    ViscousManeuveringDerivatives

Dimensional quasi-steady viscous corrections to the four linear sway-yaw
velocity derivatives. Their units and sign conventions match the
corresponding fields in [`WangHydrodynamicDerivatives`](@ref).
"""
struct ViscousManeuveringDerivatives{T}
    Y_v::T
    Y_r::T
    N_v::T
    N_r::T
end

function ViscousManeuveringDerivatives(Y_v, Y_r, N_v, N_r)
    values = promote(Y_v, Y_r, N_v, N_r)
    return ViscousManeuveringDerivatives{typeof(first(values))}(values...)
end

raw"""
    ViscousManeuveringCorrectionResult

Result of [`viscous_maneuvering_correction`](@ref). `pressure_derivatives`
come from the displacement-thickness transpiration potential and
`shear_derivatives` come directly from wall traction. Their sum is stored in
`derivatives`. `coupling_iterations=0` denotes a one-way interaction; positive
values denote under-relaxed viscous-inviscid feedback passes.
"""
struct ViscousManeuveringCorrectionResult{D,B,T,P,G,R,C}
    derivatives::D
    pressure_derivatives::D
    shear_derivatives::D
    base_boundary_layer::B
    transpiration_derivatives::T
    correction_potentials::P
    correction_gradients::G
    boundary_residuals::R
    linearization::Symbol
    coupling_iterations::Int
    coupling_residual::C
end

raw"""
    head_kinetic_shape_factor(H)

Return Head's entrainment shape parameter ``H_1(H)`` using the standard
piecewise turbulent-profile correlation.
"""
function head_kinetic_shape_factor(H::Real)
    H > 1.1 || throw(DomainError(H, "Head's shape factor H must exceed 1.1"))
    if H <= 1.6
        return 3.3 + 0.8234 * (H - 1.1)^(-1.287)
    end
    return 3.3 + 1.5501 * (H - 0.6778)^(-3.064)
end

function _head_kinetic_shape_factor_derivative(H::Real)
    H > 1.1 || throw(DomainError(H, "Head's shape factor H must exceed 1.1"))
    if H <= 1.6
        return -1.287 * 0.8234 * (H - 1.1)^(-2.287)
    end
    return -3.064 * 1.5501 * (H - 0.6778)^(-4.064)
end

raw"""
    head_skin_friction_coefficient(H, Re_theta)

Return the local skin-friction coefficient from the Ludwieg-Tillmann
correlation used with Head's method,

```math
C_f=0.246\,10^{-0.678H}Re_\theta^{-0.268}.
```
"""
function head_skin_friction_coefficient(H::Real, Re_theta::Real)
    H > 1.1 || throw(DomainError(H, "Head's shape factor H must exceed 1.1"))
    Re_theta > 0 || throw(DomainError(
        Re_theta,
        "momentum-thickness Reynolds number must be positive",
    ))
    return 0.246 * 10^(-0.678 * H) * Re_theta^(-0.268)
end

raw"""
    head_entrainment_coefficient(H1)

Return Head's turbulent entrainment coefficient,

```math
E=0.0306(H_1-3)^{-0.6169}.
```
"""
function head_entrainment_coefficient(H1::Real)
    H1 > 3 || throw(DomainError(H1, "Head's H1 must exceed 3"))
    return 0.0306 * (H1 - 3)^(-0.6169)
end

function _head_rhs(theta, H, edge_speed, edge_speed_gradient, viscosity)
    H1 = head_kinetic_shape_factor(H)
    dH1_dH = _head_kinetic_shape_factor_derivative(H)
    Re_theta = edge_speed * theta / viscosity
    skin_friction = head_skin_friction_coefficient(H, Re_theta)
    entrainment = head_entrainment_coefficient(H1)
    logarithmic_speed_gradient = edge_speed_gradient / edge_speed
    theta_gradient = skin_friction / 2 -
        (H + 2) * theta * logarithmic_speed_gradient
    shape_gradient = (
        entrainment - H1 * theta_gradient -
        theta * H1 * logarithmic_speed_gradient
    ) / (theta * dH1_dH)
    return (theta_gradient, shape_gradient)
end

function _rk4_head_step(theta, H, U_left, U_right, ds, viscosity)
    dU_ds = (U_right - U_left) / ds
    U_mid = (U_left + U_right) / 2
    k1_theta, k1_H = _head_rhs(theta, H, U_left, dU_ds, viscosity)
    k2_theta, k2_H = _head_rhs(
        theta + ds * k1_theta / 2,
        H + ds * k1_H / 2,
        U_mid,
        dU_ds,
        viscosity,
    )
    k3_theta, k3_H = _head_rhs(
        theta + ds * k2_theta / 2,
        H + ds * k2_H / 2,
        U_mid,
        dU_ds,
        viscosity,
    )
    k4_theta, k4_H = _head_rhs(
        theta + ds * k3_theta,
        H + ds * k3_H,
        U_right,
        dU_ds,
        viscosity,
    )
    next_theta = theta + ds * (
        k1_theta + 2k2_theta + 2k3_theta + k4_theta
    ) / 6
    next_H = H + ds * (k1_H + 2k2_H + 2k3_H + k4_H) / 6
    return next_theta, next_H
end

raw"""
    solve_head_boundary_layer(
        coordinate,
        edge_speed,
        kinematic_viscosity;
        closure=HeadTurbulentClosure(),
        initial_momentum_thickness=nothing,
        integration_substeps=8,
        stop_at_separation=true,
    )

March Head's two-equation turbulent integral boundary layer along a prescribed
external-speed distribution. The governing equations are

```math
\frac{d\theta}{ds}+(H+2)\frac{\theta}{U_e}\frac{dU_e}{ds}
=\frac{C_f}{2},
\qquad
\frac{d}{ds}(U_e\theta H_1)=U_eE.
```

Coordinates must be positive and strictly increasing, and edge speeds must be
positive. If no initial momentum thickness is supplied, the first station is
initialized with the turbulent flat-plate estimate
``\theta=0.036sRe_s^{-1/5}``. This formulation assumes a fully turbulent,
attached layer. `separated` is therefore a model-validity flag, not a
post-separation solution. `integration_substeps` stabilizes the march across
coarse external-speed intervals. With `stop_at_separation=true`, the state is
frozen once `H` reaches the separation threshold.
"""
function solve_head_boundary_layer(
    coordinate::AbstractVector,
    edge_speed::AbstractVector,
    kinematic_viscosity::Real;
    closure::HeadTurbulentClosure=HeadTurbulentClosure(),
    initial_momentum_thickness=nothing,
    integration_substeps::Integer=8,
    stop_at_separation::Bool=true,
)
    station_count = length(coordinate)
    station_count > 0 || throw(ArgumentError("at least one station is required"))
    length(edge_speed) == station_count || throw(DimensionMismatch(
        "coordinate and edge_speed must have equal length",
    ))
    kinematic_viscosity > 0 || throw(ArgumentError(
        "kinematic_viscosity must be positive",
    ))
    coordinate[1] > 0 || throw(ArgumentError(
        "the first boundary-layer coordinate must be positive",
    ))
    all(diff(coordinate) .> 0) || throw(ArgumentError(
        "boundary-layer coordinates must be strictly increasing",
    ))
    all(edge_speed .> 0) || throw(ArgumentError("edge speeds must be positive"))
    integration_substeps > 0 || throw(ArgumentError(
        "integration_substeps must be positive",
    ))

    initial_theta = if isnothing(initial_momentum_thickness)
        first_Reynolds_number = edge_speed[1] * coordinate[1] / kinematic_viscosity
        0.036 * coordinate[1] * first_Reynolds_number^(-1 / 5)
    else
        initial_momentum_thickness > 0 || throw(ArgumentError(
            "initial_momentum_thickness must be positive",
        ))
        initial_momentum_thickness
    end
    T = promote_type(
        eltype(coordinate),
        eltype(edge_speed),
        typeof(kinematic_viscosity),
        typeof(closure.initial_shape_factor),
        typeof(initial_theta),
    )
    theta = Vector{T}(undef, station_count)
    H = Vector{T}(undef, station_count)
    H1 = Vector{T}(undef, station_count)
    displacement = Vector{T}(undef, station_count)
    skin_friction = Vector{T}(undef, station_count)
    separated = falses(station_count)
    theta[1] = initial_theta
    H[1] = closure.initial_shape_factor
    separation_reached = H[1] >= closure.separation_shape_factor

    for station in 2:station_count
        if stop_at_separation && separation_reached
            theta[station] = theta[station - 1]
            H[station] = H[station - 1]
            continue
        end
        ds = coordinate[station] - coordinate[station - 1]
        substep_length = ds / integration_substeps
        speed_change = edge_speed[station] - edge_speed[station - 1]
        next_theta = theta[station - 1]
        next_H = H[station - 1]
        for substep in 1:integration_substeps
            left_fraction = (substep - 1) / integration_substeps
            right_fraction = substep / integration_substeps
            U_left = edge_speed[station - 1] + left_fraction * speed_change
            U_right = edge_speed[station - 1] + right_fraction * speed_change
            next_theta, next_H = _rk4_head_step(
                next_theta,
                next_H,
                U_left,
                U_right,
                substep_length,
                kinematic_viscosity,
            )
            if stop_at_separation &&
                    next_H >= closure.separation_shape_factor
                next_H = oftype(next_H, closure.separation_shape_factor)
                separation_reached = true
                break
            end
        end
        next_theta > 0 || throw(DomainError(
            next_theta,
            "Head integration produced nonpositive momentum thickness",
        ))
        next_H > 1.1 || throw(DomainError(
            next_H,
            "Head integration left the valid turbulent-profile range H > 1.1",
        ))
        theta[station] = next_theta
        H[station] = next_H
    end

    for station in eachindex(theta)
        H1[station] = head_kinetic_shape_factor(H[station])
        displacement[station] = H[station] * theta[station]
        Re_theta = edge_speed[station] * theta[station] / kinematic_viscosity
        skin_friction[station] = head_skin_friction_coefficient(
            H[station],
            Re_theta,
        )
        separated[station] = H[station] >= closure.separation_shape_factor
        if stop_at_separation && separated[station]
            skin_friction[station] = zero(skin_friction[station])
        end
    end
    return IntegralBoundaryLayerResult(
        coordinate,
        edge_speed,
        theta,
        H,
        H1,
        displacement,
        skin_friction,
        separated,
    )
end

function _strip_coordinates_and_lengths(centers::AbstractMatrix)
    station_count = size(centers, 1)
    station_count > 0 || throw(ArgumentError("a strip must contain at least one panel"))
    T = eltype(centers)
    coordinate = Vector{T}(undef, station_count)
    panel_length = Vector{T}(undef, station_count)
    if station_count == 1
        coordinate[1] = one(T)
        panel_length[1] = one(T)
        return coordinate, panel_length
    end

    spacing = Vector{T}(undef, station_count - 1)
    for station in eachindex(spacing)
        spacing[station] = norm(
            @view(centers[station + 1, :]) - @view(centers[station, :]),
        )
        spacing[station] > 0 || throw(ArgumentError(
            "successive strip-panel centers must be distinct",
        ))
    end
    coordinate[1] = spacing[1] / 2
    for station in 2:station_count
        coordinate[station] = coordinate[station - 1] + spacing[station - 1]
    end
    panel_length[1] = spacing[1]
    panel_length[end] = spacing[end]
    for station in 2:(station_count - 1)
        panel_length[station] = (spacing[station - 1] + spacing[station]) / 2
    end
    return coordinate, panel_length
end

function _mass_defect_divergence(panel_length, width, mass_defect)
    station_count = length(panel_length)
    length(width) == station_count || throw(DimensionMismatch(
        "panel_length and width must have equal length",
    ))
    length(mass_defect) == station_count || throw(DimensionMismatch(
        "panel_length and mass_defect must have equal length",
    ))
    T = promote_type(eltype(panel_length), eltype(width), eltype(mass_defect))
    transpiration = Vector{T}(undef, station_count)
    face_mass_defect = Vector{T}(undef, station_count + 1)
    face_mass_defect[1] = zero(T)
    for face in 2:station_count
        face_mass_defect[face] = (
            mass_defect[face - 1] + mass_defect[face]
        ) / 2
    end
    # The outgoing mass defect is continued into the unresolved wake. This
    # finite-volume boundary condition conserves the strip mass deficit and
    # avoids making the last-panel derivative depend on a one-sided stencil.
    face_mass_defect[end] = mass_defect[end]
    for station in 1:station_count
        transpiration[station] = (
            face_mass_defect[station + 1] - face_mass_defect[station]
        ) / (
            width[station] * panel_length[station]
        )
    end
    return transpiration
end

raw"""
    body_relative_edge_velocity(
        mesh,
        forward_speed,
        sway_velocity,
        yaw_rate;
        surge_gradient=nothing,
        sway_gradient=nothing,
        yaw_gradient=nothing,
        x_reference=0,
        y_reference=0,
    )

Assemble and tangentially project the body-relative external velocity used by
the boundary-layer solver. With unit rigid-body potential gradients supplied,

```math
\boldsymbol{u}_e=
U\nabla\phi_U+v\nabla\phi_v+r\nabla\phi_r
-\left[U\boldsymbol{e}_x+v\boldsymbol{e}_y
+r\boldsymbol{e}_z\mathbin{\boldsymbol{\times}}
(\boldsymbol{x}-\boldsymbol{x}_0)\right].
```

The result is projected by ``\boldsymbol{I}-\boldsymbol{n}\boldsymbol{n}^T``
to remove panel-quadrature residuals. Omitted potential gradients contribute
zero, which provides a kinematic surface-projection approximation.
"""
function body_relative_edge_velocity(
    mesh::Mesh,
    forward_speed::Real,
    sway_velocity::Real,
    yaw_rate::Real;
    surge_gradient::Union{Nothing,AbstractMatrix}=nothing,
    sway_gradient::Union{Nothing,AbstractMatrix}=nothing,
    yaw_gradient::Union{Nothing,AbstractMatrix}=nothing,
    x_reference::Real=0,
    y_reference::Real=0,
)
    forward_speed > 0 || throw(ArgumentError("forward_speed must be positive"))
    gradients = (
        surge_gradient=surge_gradient,
        sway_gradient=sway_gradient,
        yaw_gradient=yaw_gradient,
    )
    for (name, gradient) in pairs(gradients)
        if !isnothing(gradient)
            size(gradient) == (mesh.nfaces, 3) || throw(DimensionMismatch(
                "$(name) must have size (mesh.nfaces, 3)",
            ))
        end
    end
    gradient_types = (
        eltype(gradient) for gradient in values(gradients) if !isnothing(gradient)
    )
    T = promote_type(
        eltype(mesh.centers),
        typeof(forward_speed),
        typeof(sway_velocity),
        typeof(yaw_rate),
        typeof(x_reference),
        typeof(y_reference),
        gradient_types...,
    )
    edge_velocity = Matrix{T}(undef, mesh.nfaces, 3)
    for panel in 1:mesh.nfaces
        x = mesh.centers[panel, 1] - x_reference
        y = mesh.centers[panel, 2] - y_reference
        body_x = forward_speed - yaw_rate * y
        body_y = sway_velocity + yaw_rate * x
        potential_x = zero(T)
        potential_y = zero(T)
        potential_z = zero(T)
        if !isnothing(surge_gradient)
            potential_x += forward_speed * surge_gradient[panel, 1]
            potential_y += forward_speed * surge_gradient[panel, 2]
            potential_z += forward_speed * surge_gradient[panel, 3]
        end
        if !isnothing(sway_gradient)
            potential_x += sway_velocity * sway_gradient[panel, 1]
            potential_y += sway_velocity * sway_gradient[panel, 2]
            potential_z += sway_velocity * sway_gradient[panel, 3]
        end
        if !isnothing(yaw_gradient)
            potential_x += yaw_rate * yaw_gradient[panel, 1]
            potential_y += yaw_rate * yaw_gradient[panel, 2]
            potential_z += yaw_rate * yaw_gradient[panel, 3]
        end
        relative_x = potential_x - body_x
        relative_y = potential_y - body_y
        relative_z = potential_z
        normal_component =
            relative_x * mesh.normals[panel, 1] +
            relative_y * mesh.normals[panel, 2] +
            relative_z * mesh.normals[panel, 3]
        edge_velocity[panel, 1] =
            relative_x - normal_component * mesh.normals[panel, 1]
        edge_velocity[panel, 2] =
            relative_y - normal_component * mesh.normals[panel, 2]
        edge_velocity[panel, 3] =
            relative_z - normal_component * mesh.normals[panel, 3]
    end
    return edge_velocity
end

raw"""
    solve_quasi3d_boundary_layer(
        grid,
        edge_velocity,
        kinematic_viscosity;
        rho=SETTINGS.rho,
        reference=(0, 0, 0),
        closure=HeadTurbulentClosure(),
        minimum_edge_speed=nothing,
        integration_substeps=8,
        stop_at_separation=true,
    )

March a Head integral boundary layer along every longitudinal strip of a
[`StructuredPanelGrid`](@ref). The scalar development uses ``|\boldsymbol{u}_e|``;
the wall-shear direction retains all three components of the tangential edge
velocity. The strip-width approximation is ``w_i=\Delta S_i/\Delta s_i`` and
the displacement transpiration is

```math
v_n^{BL}=\frac{1}{w}\frac{d}{ds}
\left(wU_e\delta^*\right).
```

This is a quasi-three-dimensional attached-flow model. It does not solve the
crossflow momentum-integral equation and does not model a separated bilge
vortex or the downstream wake.
"""
function solve_quasi3d_boundary_layer(
    grid::StructuredPanelGrid,
    edge_velocity::AbstractMatrix,
    kinematic_viscosity::Real;
    rho::Real=SETTINGS.rho,
    reference=(0, 0, 0),
    closure::HeadTurbulentClosure=HeadTurbulentClosure(),
    minimum_edge_speed=nothing,
    integration_substeps::Integer=8,
    stop_at_separation::Bool=true,
)
    mesh = grid.mesh
    size(edge_velocity) == (mesh.nfaces, 3) || throw(DimensionMismatch(
        "edge_velocity must have size (mesh.nfaces, 3)",
    ))
    kinematic_viscosity > 0 || throw(ArgumentError(
        "kinematic_viscosity must be positive",
    ))
    rho > 0 || throw(ArgumentError("rho must be positive"))
    length(reference) == 3 || throw(DimensionMismatch(
        "reference must contain three coordinates",
    ))
    panel_counts = zeros(Int, mesh.nfaces)
    for strip in grid.strips, panel in strip
        1 <= panel <= mesh.nfaces || throw(BoundsError(mesh.areas, panel))
        panel_counts[panel] += 1
    end
    all(panel_counts .== 1) || throw(ArgumentError(
        "grid strips must contain every mesh panel exactly once",
    ))

    raw_edge_speed = [norm(@view(edge_velocity[panel, :])) for panel in 1:mesh.nfaces]
    maximum_edge_speed = maximum(raw_edge_speed)
    maximum_edge_speed > 0 || throw(ArgumentError(
        "at least one panel must have nonzero edge speed",
    ))
    speed_floor = isnothing(minimum_edge_speed) ?
        sqrt(eps(Float64)) * maximum_edge_speed : minimum_edge_speed
    speed_floor > 0 || throw(ArgumentError("minimum_edge_speed must be positive"))
    T = promote_type(
        eltype(edge_velocity),
        eltype(mesh.centers),
        eltype(mesh.areas),
        typeof(kinematic_viscosity),
        typeof(rho),
        typeof(speed_floor),
    )
    edge_speed = Vector{T}(undef, mesh.nfaces)
    momentum_thickness = Vector{T}(undef, mesh.nfaces)
    displacement_thickness = Vector{T}(undef, mesh.nfaces)
    shape_factor = Vector{T}(undef, mesh.nfaces)
    skin_friction = Vector{T}(undef, mesh.nfaces)
    wall_shear = Matrix{T}(undef, mesh.nfaces, 3)
    transpiration = Vector{T}(undef, mesh.nfaces)
    separated = falses(mesh.nfaces)
    for strip in grid.strips
        centers = mesh.centers[strip, :]
        coordinate, panel_length = _strip_coordinates_and_lengths(centers)
        strip_speed = Vector{T}(undef, length(strip))
        for (local_index, panel) in enumerate(strip)
            speed = raw_edge_speed[panel]
            strip_speed[local_index] = speed > speed_floor ? speed : speed_floor
        end
        strip_result = solve_head_boundary_layer(
            coordinate,
            strip_speed,
            kinematic_viscosity;
            closure,
            integration_substeps,
            stop_at_separation,
        )
        width = mesh.areas[strip] ./ panel_length
        first_separated = findfirst(strip_result.separated)
        attached_count = isnothing(first_separated) ?
            length(strip) : first_separated - 1
        strip_transpiration = zeros(T, length(strip))
        if attached_count > 0
            attached_indices = 1:attached_count
            mass_defect = width[attached_indices] .*
                strip_result.edge_speed[attached_indices] .*
                strip_result.displacement_thickness[attached_indices]
            strip_transpiration[attached_indices] .= _mass_defect_divergence(
                panel_length[attached_indices],
                width[attached_indices],
                mass_defect,
            )
        end

        for (local_index, panel) in enumerate(strip)
            speed = strip_result.edge_speed[local_index]
            edge_speed[panel] = speed
            momentum_thickness[panel] = strip_result.momentum_thickness[local_index]
            displacement_thickness[panel] =
                strip_result.displacement_thickness[local_index]
            shape_factor[panel] = strip_result.shape_factor[local_index]
            skin_friction[panel] =
                strip_result.skin_friction_coefficient[local_index]
            transpiration[panel] = strip_transpiration[local_index]
            separated[panel] = strip_result.separated[local_index]
            shear_scale = rho * skin_friction[panel] * speed / 2
            wall_shear[panel, 1] = shear_scale * edge_velocity[panel, 1]
            wall_shear[panel, 2] = shear_scale * edge_velocity[panel, 2]
            wall_shear[panel, 3] = shear_scale * edge_velocity[panel, 3]
        end
    end

    force = zeros(T, 3)
    moment = zeros(T, 3)
    for panel in 1:mesh.nfaces
        panel_force = @view(wall_shear[panel, :]) .* mesh.areas[panel]
        force .+= panel_force
        arm = @view(mesh.centers[panel, :]) .- reference
        moment .+= cross(arm, panel_force)
    end
    return Quasi3DBoundaryLayerResult(
        edge_speed,
        momentum_thickness,
        displacement_thickness,
        shape_factor,
        skin_friction,
        wall_shear,
        transpiration,
        separated,
        force,
        moment,
    )
end

raw"""
    separation_stern_mask(grid, layer; mode=:per_strip, margin=0)

Build a velocity-integral truncation mask from a computed boundary-layer
solution, as a physically derived alternative to the semi-empirical Schmitz cut
of [`wang_stern_mask`](@ref).

Wang et al. justify their stern truncation by the vorticity shed from the after
body, which unloads the stern pressures. That argument identifies *where* the
cut belongs — at the onset of separation — but the maximum-area section is only
a geometric proxy for it. Given a [`Quasi3DBoundaryLayerResult`](@ref) or a
[`ThreeDimensionalBoundaryLayerResult`](@ref), the separation location is
predicted rather than assumed.

`mode` selects how the separation flags become a mask:

- `:per_strip` keeps, on each strip independently, the panels upstream of that
  strip's first separated station. This follows the three-dimensional
  separation line, which on a full-form hull reaches much further forward near
  the keel than near the waterline.
- `:first_station` finds the most upstream separated panel anywhere on the hull
  and applies that single longitudinal station as a plane cut, which is the
  closest like-for-like comparison with `wang_stern_mask`.
- `:none` keeps every panel, for a whole-hull reference.

`margin` shifts the cut by a number of stations: positive values move it
forward (more conservative, discarding more of the stern).

Strips whose layer never separates are kept in full.
"""
function separation_stern_mask(
    grid::StructuredPanelGrid,
    separated::AbstractVector{Bool};
    mode::Symbol=:per_strip,
    margin::Integer=0,
)
    mesh = grid.mesh
    length(separated) == mesh.nfaces || throw(DimensionMismatch(
        "separated must contain one flag per mesh face",
    ))
    mode in (:per_strip, :first_station, :none) || throw(ArgumentError(
        "mode must be :per_strip, :first_station, or :none",
    ))
    mode === :none && return trues(mesh.nfaces)

    mask = trues(mesh.nfaces)
    if mode === :per_strip
        for strip in grid.strips
            flags = [separated[panel] for panel in strip]
            onset = findfirst(flags)
            isnothing(onset) && continue
            cut = max(1, onset - margin)
            for local_index in cut:length(strip)
                mask[strip[local_index]] = false
            end
        end
        return mask
    end

    # :first_station — the most upstream separation anywhere sets one plane cut.
    earliest = nothing
    for strip in grid.strips
        flags = [separated[panel] for panel in strip]
        onset = findfirst(flags)
        isnothing(onset) && continue
        earliest = isnothing(earliest) ? onset : min(earliest, onset)
    end
    isnothing(earliest) && return mask
    cut = max(1, earliest - margin)
    cut_x = maximum(
        mesh.centers[strip[min(cut, length(strip))], 1] for strip in grid.strips
    )
    for panel in 1:mesh.nfaces
        mask[panel] = mesh.centers[panel, 1] >= cut_x
    end
    return mask
end

separation_stern_mask(grid::StructuredPanelGrid, layer; kwargs...) =
    separation_stern_mask(grid, collect(Bool, layer.separated); kwargs...)

function _velocity_derivative_sum(first::ViscousManeuveringDerivatives, second)
    return ViscousManeuveringDerivatives(
        first.Y_v + second.Y_v,
        first.Y_r + second.Y_r,
        first.N_v + second.N_v,
        first.N_r + second.N_r,
    )
end

function _tangential_projection(mesh::Mesh, velocity::AbstractMatrix)
    size(velocity) == (mesh.nfaces, 3) || throw(DimensionMismatch(
        "velocity must have size (mesh.nfaces, 3)",
    ))
    projected = similar(velocity)
    for panel in 1:mesh.nfaces
        normal_component =
            velocity[panel, 1] * mesh.normals[panel, 1] +
            velocity[panel, 2] * mesh.normals[panel, 2] +
            velocity[panel, 3] * mesh.normals[panel, 3]
        for axis in 1:3
            projected[panel, axis] = velocity[panel, axis] -
                normal_component * mesh.normals[panel, axis]
        end
    end
    return projected
end

raw"""
    _layer_geometry(closure, mesh)

Geometry a boundary-layer closure needs rebuilt once rather than per solve.

The quasi-3D march works off the grid's strips and needs nothing. The surface
solver needs the welded topology and its metrics, and a correction evaluates the
layer `coupling_iterations + 1` times for each of `1 + 2n` states under central
differencing, so rebuilding the topology inside each of those would dominate the
cost. The topology is integer connectivity and the metrics depend only on the
mesh, neither of which the sway velocity or the yaw rate touch.
"""
_layer_geometry(::HeadTurbulentClosure, ::Mesh) = nothing

function _layer_geometry(::ThreeDimensionalClosure, mesh::Mesh)
    topology = build_surface_topology(mesh)
    return (topology = topology, metrics = build_surface_metrics(mesh, topology))
end

"""
    _solve_layer(closure, grid, edge_velocity, kinematic_viscosity; kwargs...)

Solve the boundary layer with whichever method `closure` selects, presenting the
one interface the coupling needs.

Both methods return a result carrying `transpiration_velocity`, `force`,
`moment` and `separated`, which is everything
[`viscous_maneuvering_correction`](@ref) consumes.
"""
function _solve_layer(closure::HeadTurbulentClosure, grid, edge_velocity,
        kinematic_viscosity; rho, reference, minimum_edge_speed,
        integration_substeps, stop_at_separation, geometry = nothing,
        initial = nothing)
    return solve_quasi3d_boundary_layer(grid, edge_velocity, kinematic_viscosity;
        rho, reference, closure, minimum_edge_speed, integration_substeps,
        stop_at_separation)
end

function _solve_layer(closure::ThreeDimensionalClosure, grid, edge_velocity,
        kinematic_viscosity; rho, reference, minimum_edge_speed,
        integration_substeps = nothing, stop_at_separation = true,
        geometry = nothing, initial = nothing)
    topology = isnothing(geometry) ? nothing : geometry.topology
    metrics = isnothing(geometry) ? nothing : geometry.metrics
    return solve_surface_boundary_layer(grid, edge_velocity, kinematic_viscosity;
        rho, reference, closure, minimum_edge_speed, topology, metrics, initial)
end

# The converged state of a surface solve, for warm-starting the next one. The
# march has nothing to carry over.
_layer_warm_start(::HeadTurbulentClosure, result) = nothing
_layer_warm_start(::ThreeDimensionalClosure, result) = result.diagnostics.states

raw"""
    viscous_maneuvering_correction(
        grid,
        surge_gradient,
        sway_gradient,
        yaw_gradient,
        forward_speed,
        kinematic_viscosity;
        rho=SETTINGS.rho,
        x_reference=0,
        y_reference=0,
        green_functions=(Rankine(), RankineReflected()),
        integration_substeps=8,
        stop_at_separation=true,
        linearization=:auto,
        closure=HeadTurbulentClosure(),
        coupling_iterations=0,
        coupling_relaxation=0.5,
        coupling_tolerance=1e-6,
    )

Linearize the integral boundary layer with respect to sway velocity and yaw
rate. Direct shear-load derivatives are combined with a displacement-thickness
correction.

`closure` selects the method. `HeadTurbulentClosure` — the default, and what
every existing result was produced with — runs the quasi-3D strip march.
`ThreeDimensionalClosure` runs the globally coupled surface solve of
[`solve_surface_boundary_layer`](@ref), which carries a crossflow: its
transpiration is the surface divergence of ``u_e\boldsymbol{\delta}^*`` with
*both* components of the displacement thickness, where the march has only the
streamwise one, and its wall shear is skewed out of the local stream.
For each motion ``j``, the differentiated transpiration condition is solved by
the existing BEM and contributes

```math
\Delta Q_j^{p}=-\rho U\int_{S_h}
\frac{\partial\psi_j}{\partial x}g_Q\,dS.
```

This correction is intended to augment a whole-hull inviscid result. Applying
it after a Schmitz stern truncation mixes two different viscous corrections
and is not recommended. The default `coupling_iterations=0` is one-way. Each
positive coupling iteration solves the transpiration potential, projects its
surface velocity tangentially, and feeds it back to the boundary-layer edge
condition with `coupling_relaxation`, stopping once the relative change falls
below `coupling_tolerance`. With `ThreeDimensionalClosure` each iteration warm
starts from the previous converged state, which is what keeps a coupled sweep
affordable: without it the cost is `coupling_iterations * (1 + 2 * n_states)`
cold solves.
"""
function viscous_maneuvering_correction(
    grid::StructuredPanelGrid,
    surge_gradient::AbstractMatrix,
    sway_gradient::AbstractMatrix,
    yaw_gradient::AbstractMatrix,
    forward_speed::Real,
    kinematic_viscosity::Real;
    rho::Real=SETTINGS.rho,
    x_reference::Real=0,
    y_reference::Real=0,
    closure::Union{HeadTurbulentClosure,ThreeDimensionalClosure}=HeadTurbulentClosure(),
    minimum_edge_speed=nothing,
    integration_substeps::Integer=8,
    stop_at_separation::Bool=true,
    linearization::Symbol=:auto,
    relative_difference_step::Real=1e-5,
    coupling_iterations::Integer=0,
    coupling_relaxation::Real=0.5,
    coupling_tolerance::Real=1e-6,
    green_functions=(Rankine(), RankineReflected()),
)
    mesh = grid.mesh
    forward_speed > 0 || throw(ArgumentError("forward_speed must be positive"))
    kinematic_viscosity > 0 || throw(ArgumentError(
        "kinematic_viscosity must be positive",
    ))
    rho > 0 || throw(ArgumentError("rho must be positive"))
    linearization in (:auto, :forwarddiff, :central) || throw(ArgumentError(
        "linearization must be :auto, :forwarddiff, or :central",
    ))
    relative_difference_step > 0 || throw(ArgumentError(
        "relative_difference_step must be positive",
    ))
    coupling_iterations >= 0 || throw(ArgumentError(
        "coupling_iterations must be nonnegative",
    ))
    0 < coupling_relaxation <= 1 || throw(ArgumentError(
        "coupling_relaxation must be in (0, 1]",
    ))
    coupling_tolerance > 0 || throw(ArgumentError(
        "coupling_tolerance must be positive",
    ))

    wavenumber = zero(mesh.centers[1, 1])
    S, K = assemble_matrices(
        green_functions,
        mesh,
        wavenumber;
        direct=false,
    )
    layer_geometry = _layer_geometry(closure, mesh)

    function coupled_boundary_layer_state(state)
        rigid_edge_velocity = body_relative_edge_velocity(
            mesh,
            forward_speed,
            state[1],
            state[2];
            surge_gradient,
            sway_gradient,
            yaw_gradient,
            x_reference,
            y_reference,
        )
        correction_edge_velocity = zeros(
            eltype(rigid_edge_velocity),
            mesh.nfaces,
            3,
        )
        iteration_residual = zero(eltype(rigid_edge_velocity))
        warm_start = nothing
        for _ in 1:coupling_iterations
            edge_velocity = rigid_edge_velocity .+ correction_edge_velocity
            result = _solve_layer(
                closure,
                grid,
                edge_velocity,
                kinematic_viscosity;
                rho,
                reference=(x_reference, y_reference, zero(x_reference)),
                minimum_edge_speed,
                integration_substeps,
                stop_at_separation,
                geometry=layer_geometry,
                initial=warm_start,
            )
            warm_start = _layer_warm_start(closure, result)
            _, correction_sources = solve(
                K,
                S,
                result.transpiration_velocity;
                direct=false,
            )
            correction_gradient = evaluate_indirect_potential_gradient(
                green_functions,
                mesh,
                correction_sources,
                wavenumber,
            )
            proposed_edge_velocity = _tangential_projection(
                mesh,
                real.(correction_gradient),
            )
            relaxed_edge_velocity =
                (1 - coupling_relaxation) .* correction_edge_velocity .+
                coupling_relaxation .* proposed_edge_velocity
            denominator = norm(relaxed_edge_velocity)
            change = norm(relaxed_edge_velocity .- correction_edge_velocity)
            iteration_residual = iszero(denominator) ? change :
                change / denominator
            correction_edge_velocity = relaxed_edge_velocity
            # The loop computed this every pass and never acted on it, so a
            # converged coupling paid for every remaining iteration.
            iteration_residual < coupling_tolerance && break
        end
        final_edge_velocity = rigid_edge_velocity .+ correction_edge_velocity
        result = _solve_layer(
            closure,
            grid,
            final_edge_velocity,
            kinematic_viscosity;
            rho,
            reference=(x_reference, y_reference, zero(x_reference)),
            minimum_edge_speed,
            integration_substeps,
            stop_at_separation,
            geometry=layer_geometry,
            initial=warm_start,
        )
        return result, iteration_residual
    end

    function boundary_layer_outputs(state)
        result, _ = coupled_boundary_layer_state(state)
        return vcat(
            result.transpiration_velocity,
            result.force[2],
            result.moment[3],
        )
    end

    base_state = [zero(forward_speed), zero(forward_speed)]
    base_boundary_layer, base_coupling_residual =
        coupled_boundary_layer_state(base_state)
    selected_linearization = if linearization === :auto
        any(base_boundary_layer.separated) || coupling_iterations > 0 ?
            :central : :forwarddiff
    else
        linearization
    end
    if selected_linearization === :forwarddiff && coupling_iterations > 0
        throw(ArgumentError(
            "forward-mode linearization through coupled complex BEM solves " *
            "is not supported; use linearization=:central or :auto",
        ))
    end
    output_jacobian = if selected_linearization === :forwarddiff
        ForwardDiff.jacobian(boundary_layer_outputs, base_state)
    else
        output_at_base = boundary_layer_outputs(base_state)
        jacobian = Matrix{eltype(output_at_base)}(
            undef,
            length(output_at_base),
            length(base_state),
        )
        length_scale = maximum(mesh.centers[:, 1]) - minimum(mesh.centers[:, 1])
        length_scale > 0 || throw(ArgumentError(
            "mesh must have nonzero longitudinal extent",
        ))
        steps = (
            relative_difference_step * forward_speed,
            relative_difference_step * forward_speed / length_scale,
        )
        for variable in eachindex(base_state)
            positive_state = copy(base_state)
            negative_state = copy(base_state)
            positive_state[variable] += steps[variable]
            negative_state[variable] -= steps[variable]
            jacobian[:, variable] .= (
                boundary_layer_outputs(positive_state) .-
                boundary_layer_outputs(negative_state)
            ) ./ (2 * steps[variable])
        end
        jacobian
    end

    transpiration_derivatives = output_jacobian[1:mesh.nfaces, :]
    correction_potentials, correction_sources = solve(
        K,
        S,
        transpiration_derivatives;
        direct=false,
    )
    sway_correction_gradient = evaluate_indirect_potential_gradient(
        green_functions,
        mesh,
        correction_sources[:, 1],
        wavenumber,
    )
    yaw_correction_gradient = evaluate_indirect_potential_gradient(
        green_functions,
        mesh,
        correction_sources[:, 2],
        wavenumber,
    )
    correction_gradients = (
        sway=real.(sway_correction_gradient),
        yaw=real.(yaw_correction_gradient),
    )
    real_correction_potentials = real.(correction_potentials)

    # Wang et al. truncate the velocity-derivative integrals where the flow
    # leaves the after body, and the same argument applies with more force
    # here: downstream of separation the attached-flow closure is invalid, so
    # any pressure it predicts there is meaningless. The correction is
    # therefore integrated over the attached region only. The shear
    # contribution is already attached-only, because the closure sets the wall
    # stress to zero once a station separates.
    attached = .!base_boundary_layer.separated
    sway_weight = mesh.normals[:, 2] .* mesh.areas .* attached
    yaw_mode = maneuvering_boundary_condition(
        mesh,
        :yaw;
        x_reference,
        y_reference,
    )
    yaw_weight = yaw_mode .* mesh.areas .* attached
    pressure_derivatives = ViscousManeuveringDerivatives(
        -rho * forward_speed * dot(correction_gradients.sway[:, 1], sway_weight),
        -rho * forward_speed * dot(correction_gradients.yaw[:, 1], sway_weight),
        -rho * forward_speed * dot(correction_gradients.sway[:, 1], yaw_weight),
        -rho * forward_speed * dot(correction_gradients.yaw[:, 1], yaw_weight),
    )
    shear_derivatives = ViscousManeuveringDerivatives(
        output_jacobian[mesh.nfaces + 1, 1],
        output_jacobian[mesh.nfaces + 1, 2],
        output_jacobian[mesh.nfaces + 2, 1],
        output_jacobian[mesh.nfaces + 2, 2],
    )
    derivatives = _velocity_derivative_sum(
        pressure_derivatives,
        shear_derivatives,
    )

    function relative_residual(column)
        denominator = norm(transpiration_derivatives[:, column])
        residual = norm(
            K * correction_sources[:, column] -
            transpiration_derivatives[:, column],
        )
        return iszero(denominator) ? residual : residual / denominator
    end
    return ViscousManeuveringCorrectionResult(
        derivatives,
        pressure_derivatives,
        shear_derivatives,
        base_boundary_layer,
        transpiration_derivatives,
        real_correction_potentials,
        correction_gradients,
        (sway=relative_residual(1), yaw=relative_residual(2)),
        selected_linearization,
        coupling_iterations,
        base_coupling_residual,
    )
end

function viscous_maneuvering_correction(
    grid::StructuredPanelGrid,
    surge::RigidBodyPotentialResult,
    maneuvering::PotentialFlowManeuveringResult,
    forward_speed::Real,
    kinematic_viscosity::Real;
    kwargs...,
)
    return viscous_maneuvering_correction(
        grid,
        surge.potential_gradient,
        maneuvering.potential_gradients.sway,
        maneuvering.potential_gradients.yaw,
        forward_speed,
        kinematic_viscosity;
        kwargs...,
    )
end

raw"""
    apply_viscous_correction(inviscid, correction)

Add quasi-steady viscous velocity derivatives to an inviscid
[`WangHydrodynamicDerivatives`](@ref), leaving the acceleration derivatives
unchanged.
"""
function apply_viscous_correction(
    inviscid::WangHydrodynamicDerivatives,
    correction::ViscousManeuveringDerivatives,
)
    return WangHydrodynamicDerivatives(
        inviscid.Y_vdot,
        inviscid.Y_v + correction.Y_v,
        inviscid.Y_rdot,
        inviscid.Y_r + correction.Y_r,
        inviscid.N_vdot,
        inviscid.N_v + correction.N_v,
        inviscid.N_rdot,
        inviscid.N_r + correction.N_r,
    )
end

apply_viscous_correction(
    inviscid::WangHydrodynamicDerivatives,
    correction::ViscousManeuveringCorrectionResult,
) = apply_viscous_correction(inviscid, correction.derivatives)

raw"""
    nondimensionalize_viscous_derivatives(
        derivatives,
        length,
        forward_speed;
        rho=SETTINGS.rho,
    )

Apply Wang's velocity-derivative scales to a
[`ViscousManeuveringDerivatives`](@ref).
"""
function nondimensionalize_viscous_derivatives(
    derivatives::ViscousManeuveringDerivatives,
    length::Real,
    forward_speed::Real;
    rho::Real=SETTINGS.rho,
)
    length > 0 || throw(ArgumentError("length must be positive"))
    forward_speed > 0 || throw(ArgumentError("forward_speed must be positive"))
    rho > 0 || throw(ArgumentError("rho must be positive"))
    half_rho = rho / 2
    return ViscousManeuveringDerivatives(
        derivatives.Y_v / (half_rho * length^2 * forward_speed),
        derivatives.Y_r / (half_rho * length^3 * forward_speed),
        derivatives.N_v / (half_rho * length^3 * forward_speed),
        derivatives.N_r / (half_rho * length^4 * forward_speed),
    )
end
