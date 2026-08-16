raw"""
    SectionalCrossflowDrag

Per-station geometry and separation state used by the cross-flow drag model.

`station_x` runs along the hull, `draft` is each section's immersed depth,
`spacing` the length each station stands for, and `separated_fraction` the share
of that section's panels the boundary layer reports as separated.
"""
struct SectionalCrossflowDrag{T}
    station_x::Vector{T}
    draft::Vector{T}
    spacing::Vector{T}
    separated_fraction::Vector{T}
end

raw"""
    sectional_crossflow_geometry(grid, separated)

Reduce a structured hull and a per-panel separation flag to the sectional
quantities the cross-flow drag model needs.

A strip runs longitudinally, so a station is one index into every strip at once.
The draft is taken from the panel centres of the station, which understates the
true section depth by half a panel; that bias is common to every station and so
falls out of the *shape* of the loading, not its distribution.
"""
function sectional_crossflow_geometry(grid::StructuredPanelGrid,
        separated::AbstractVector{Bool})
    mesh = grid.mesh
    isempty(grid.strips) && throw(ArgumentError("grid has no strips"))
    stations = minimum(length, grid.strips)
    stations > 1 || throw(ArgumentError("grid needs at least two stations"))

    element_type = eltype(mesh.centers)
    station_x = zeros(element_type, stations)
    draft = zeros(element_type, stations)
    spacing = zeros(element_type, stations)
    fraction = zeros(element_type, stations)
    for station in 1:stations
        panels = [strip[station] for strip in grid.strips if length(strip) >= station]
        station_x[station] = sum(mesh.centers[panel, 1] for panel in panels) /
                             length(panels)
        draft[station] = -minimum(mesh.centers[panel, 3] for panel in panels)
        fraction[station] = count(panel -> separated[panel], panels) / length(panels)
    end

    order = sortperm(station_x)
    station_x = station_x[order]
    draft = draft[order]
    fraction = fraction[order]
    for station in 1:stations
        lower = station == 1 ? station_x[1] :
                (station_x[station] + station_x[station - 1]) / 2
        upper = station == stations ? station_x[stations] :
                (station_x[station] + station_x[station + 1]) / 2
        spacing[station] = upper - lower
    end
    return SectionalCrossflowDrag(station_x, draft, spacing, fraction)
end

raw"""
    crossflow_drag_loads(
        sections,
        sway_velocity,
        yaw_rate;
        rho = SETTINGS.rho,
        drag_coefficient = 0.8,
        x_reference = 0,
    )

Sway force and yaw moment from cross-flow drag on the separated part of the
hull.

Each section is treated as a bluff two-dimensional body in the local lateral
stream, following the Fedyaevsky–Sobolev picture used throughout maneuvering
practice:

```math
Y=-\frac{\rho}{2}\int C_D(x)\,T(x)\,\left|U_c(x)\right|U_c(x)\,dx,
\qquad U_c(x)=v+r\left(x-x_0\right),
```

with the moment taking the same integrand weighted by ``x-x_0``.

What the boundary layer supplies is *where* this acts. The sectional drag
coefficient is scaled by the fraction of each section the layer reports as
separated, so a section whose flow is entirely attached contributes nothing and
the loading grows over the after body as the separation line moves forward with
drift. That is the intended replacement for a fitted stern cut: the extent comes
from the solved layer rather than from a constant chosen to match data.

`drag_coefficient` is the *published* sectional value for ship-like sections,
not a fitted one; 0.8 sits in the middle of the usual 0.5 to 1.0 range. It is a
keyword so the sensitivity to it can be reported rather than hidden.

Note what this can and cannot do. The integrand is quadratic in the lateral
velocity, so its derivative at zero drift vanishes identically: cross-flow drag
contributes to ``Y_{v|v|}`` and ``N_{v|v|}`` and contributes **nothing** to the
linear ``Y_v`` or ``N_v``. Closing the linear deficit needs the other half of the
picture — a slender-body momentum balance terminated at the separation station
rather than at the tail — which is a separate mechanism and is not this function.
"""
function crossflow_drag_loads(sections::SectionalCrossflowDrag, sway_velocity::Real,
        yaw_rate::Real; rho::Real = SETTINGS.rho, drag_coefficient::Real = 0.8,
        x_reference::Real = 0)
    rho > 0 || throw(ArgumentError("rho must be positive"))
    drag_coefficient >= 0 ||
        throw(ArgumentError("drag_coefficient must be nonnegative"))
    force = zero(promote_type(typeof(sway_velocity), typeof(yaw_rate),
        eltype(sections.draft)))
    moment = zero(force)
    for station in eachindex(sections.station_x)
        arm = sections.station_x[station] - x_reference
        lateral = sway_velocity + yaw_rate * arm
        strength = drag_coefficient * sections.separated_fraction[station] *
                   sections.draft[station] * sections.spacing[station]
        contribution = -rho / 2 * strength * abs(lateral) * lateral
        force += contribution
        moment += contribution * arm
    end
    return (force = force, moment = moment)
end

function crossflow_drag_loads(grid::StructuredPanelGrid,
        separated::AbstractVector{Bool}, args...; kwargs...)
    return crossflow_drag_loads(sectional_crossflow_geometry(grid, separated),
        args...; kwargs...)
end
