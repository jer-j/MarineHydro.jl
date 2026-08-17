raw"""
    sectional_added_mass(sections; rho = SETTINGS.rho, coefficient = 1.0)

Sway added mass per unit length of each station, ``m_{22}(x)``.

The double body of a section of draft ``T`` is a shape of depth ``2T``, whose
sway added mass is ``\rho\pi T^2`` for the flat-plate limit. `coefficient`
carries the section-shape correction — a Lewis-form factor — and is one by
default, which is the flat-plate value and close to right for the full,
box-like sections of a tanker.
"""
function sectional_added_mass(sections::SectionalCrossflowDrag;
        rho::Real = SETTINGS.rho, coefficient::Real = 1.0)
    rho > 0 || throw(ArgumentError("rho must be positive"))
    coefficient >= 0 || throw(ArgumentError("coefficient must be nonnegative"))
    return rho * π * coefficient .* sections.draft .^ 2
end

raw"""
    shed_vorticity_derivatives(
        sections, separation_x, forward_speed;
        rho = SETTINGS.rho, x_reference = 0, coefficient = 1.0,
    )

Linear sway force and yaw moment derivatives from vorticity shed at a
separation station, by the slender-body momentum balance.

A fluid section sweeps from bow to stern carrying lateral momentum
``m_{22}(x)\,v``. In steady flow the force on the hull is the rate at which that
momentum changes,

```math
f_y(x)=U\,v\,\frac{\mathrm{d}m_{22}}{\mathrm{d}x},
```

so on a *closed* body the integral telescopes to zero — d'Alembert's paradox,
and exactly the ``Y_v=0`` the potential-flow solve returns. What breaks it is
the flow leaving the surface. Downstream of the separation station ``x_s`` the
fluid is no longer bound to the hull, so the momentum ``m_{22}(x_s)v`` it holds
there is carried away into the trailing vortices instead of being recovered:

```math
Y_v=-U\,m_{22}(x_s),\qquad
N_v=-U\left[(x_s-x_0)\,m_{22}(x_s)+\int_{x_s}^{\text{bow}}m_{22}\,\mathrm{d}x\right].
```

This is the linear counterpart of [`crossflow_drag_loads`](@ref), which is
quadratic in the drift and so contributes nothing here. Between them they are
the two halves of the bilge-vortex contribution that Tanaka treats as
superposable on the boundary layer.

Two properties are worth being explicit about, because they are easy to state
wrongly.

The balance telescopes to zero — d'Alembert — only when ``m_{22}`` vanishes at
the tail, and ``m_{22}`` goes as the square of the **draft**, not the beam. A
hull whose sections taper away in depth, such as a surface-piercing ellipsoid,
recovers ``Y_v=0``. A hull carrying finite draft to its stern does not, and
slender-body theory gives it a side force whether or not the flow separates
early: the momentum simply leaves at the tail instead.

Consequently ``Y_v\simeq-2\pi(T/L)^2`` in nondimensional form, and on a hull of
near-constant after-body draft it is largely *insensitive* to where exactly the
separation station falls. That is a feature rather than a limitation — the
linear derivative does not hinge on pinning the separation line precisely — but
it also means the boundary layer's role here is to establish that the flow does
leave the after body, not to locate it to the panel.

What this does *not* do is predict ``x_s``; that comes from the boundary layer.
The whole point of the exercise is that the separation station is solved for
rather than fitted, which is what distinguishes it from a Schmitz stern cut.
"""
function shed_vorticity_derivatives(sections::SectionalCrossflowDrag,
        separation_x::Real, forward_speed::Real; rho::Real = SETTINGS.rho,
        x_reference::Real = 0, coefficient::Real = 1.0)
    forward_speed > 0 || throw(ArgumentError("forward_speed must be positive"))
    added_mass = sectional_added_mass(sections; rho, coefficient)
    stations = sections.station_x
    bow = maximum(stations)
    stern = minimum(stations)
    stern <= separation_x <= bow || throw(ArgumentError(
        "separation_x must lie between $(stern) and $(bow)"))

    # Sectional added mass at the separation station, and its integral from
    # there to the bow. Linear interpolation, so that moving the separation
    # station moves the answer smoothly — the derivative with respect to the
    # mesh has to pass through here.
    at_separation = _interpolate(stations, added_mass, separation_x)
    integral = zero(eltype(added_mass))
    for index in eachindex(stations)
        stations[index] >= separation_x || continue
        lower = max(stations[index] - sections.spacing[index] / 2, separation_x)
        upper = stations[index] + sections.spacing[index] / 2
        upper > lower || continue
        integral += added_mass[index] * (upper - lower)
    end

    sway = -forward_speed * at_separation
    yaw = -forward_speed *
          ((separation_x - x_reference) * at_separation + integral)
    return (Y_v = sway, N_v = yaw, sectional_added_mass = at_separation,
        added_mass_forward = integral)
end

# Piecewise linear lookup, clamped at the ends.
function _interpolate(stations, values, location)
    location <= first(stations) && return first(values)
    location >= last(stations) && return last(values)
    upper = findfirst(>=(location), stations)
    upper == 1 && return first(values)
    lower = upper - 1
    span = stations[upper] - stations[lower]
    span > 0 || return values[lower]
    weight = (location - stations[lower]) / span
    return values[lower] + weight * (values[upper] - values[lower])
end

raw"""
    separation_station(grid, active; quantile = 0.5)

Longitudinal station at which the boundary layer stops covering the hull.

Takes the station where the retained fraction of a section's panels first falls
below `quantile`, sweeping from the bow aft. That is the station at which the
solved layer ceases to describe the flow — because the thin-layer assumption
failed there, because the flow separated, or both — and so the station at which
its momentum is handed to the trailing vortices.
"""
function separation_station(grid::StructuredPanelGrid, active::AbstractVector{Bool};
        quantile::Real = 0.5)
    0 < quantile <= 1 || throw(ArgumentError("quantile must lie in (0, 1]"))
    mesh = grid.mesh
    sections = sectional_crossflow_geometry(grid, falses(mesh.nfaces))
    stations = sections.station_x
    retained = zeros(eltype(stations), length(stations))
    order = sortperm([sum(mesh.centers[strip, 1]) / length(strip)
                      for strip in grid.strips])
    count_at = zeros(Int, length(stations))
    kept_at = zeros(Int, length(stations))
    for strip in grid.strips, (position, panel) in enumerate(strip)
        position <= length(stations) || continue
        count_at[position] += 1
        active[panel] && (kept_at[position] += 1)
    end
    # `sectional_crossflow_geometry` sorts its stations by x; the strip index
    # runs the other way when the importer numbers from the bow, so match them
    # on the station coordinate rather than assuming an order.
    strip_x = [sum(mesh.centers[[strip[position] for strip in grid.strips
                                 if length(strip) >= position], 1]) /
               max(count_at[position], 1) for position in eachindex(stations)]
    permutation = sortperm(strip_x)
    for (slot, position) in enumerate(permutation)
        retained[slot] = count_at[position] == 0 ? one(eltype(stations)) :
                         kept_at[position] / count_at[position]
    end
    isempty(order) && return first(stations)

    for slot in length(stations):-1:1
        retained[slot] >= quantile && return stations[slot]
    end
    return first(stations)
end
