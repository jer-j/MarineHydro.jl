abstract type GreensFunction end

const _GreensFunctionVector = AbstractVector{<:GreensFunction}

function _tuple_sum(f, values::Tuple{T}) where {T}
    return f(first(values))
end

function _tuple_sum(f, values::Tuple{T,U,Vararg{Any}}) where {T,U}
    return f(first(values)) + _tuple_sum(f, Base.tail(values))
end

function _tuple_broadcast_sum(f, values::Tuple{T}) where {T}
    return f(first(values))
end

function _tuple_broadcast_sum(f, values::Tuple{T,U,Vararg{Any}}) where {T,U}
    return f(first(values)) .+ _tuple_broadcast_sum(f, Base.tail(values))
end


"""
    greens(::GreensFunction, element_1, element_2, wavenumber=nothing)

Calculates the Green's function between two points.

# Arguments
- `element_1`: first point, as a `NamedTuple` with the field `center`.
- `element_2`: second point, as a `NamedTuple` with the field `center`.
- `wavenumber`: Defines the wavenumber for waves ~~~~~. Optional for some Greens functions such as `Rankine()`.
"""
function greens end

greens(gfs::Tuple, e1, e2, w) =
    _tuple_sum(gf -> greens(gf, e1, e2, w), gfs)
greens(gfs::_GreensFunctionVector, e1, e2, w) =
    sum(greens(gf, e1, e2, w) for gf in gfs)


"""
    gradient_greens(::GreensFunction, element_1, element_2, wavenumber; with_respect_to_first_variable=false)

Calculates the gradient of the Rankine Green's function between two points.

# Arguments
- `element_1`: first point, as a `NamedTuple` with the field `center`.
- `element_2`: second point, as a `NamedTuple` with the field `center`.
- `wavenumber`: Defines the wavenumber for waves ~~~~~. Optional for some Greens functions such as `Rankine()`.
- `with_respect_to_first_variable`: A boolean flag (default is `false`). If `true`, computes the gradient with respect to `element_1`; otherwise, computes the gradient with respect to `element_2`.
"""
function gradient_greens end

gradient_greens(gfs::Tuple, e1, e2, w; with_respect_to_first_variable=false) =
    _tuple_sum(
        gf -> gradient_greens(
            gf,
            e1,
            e2,
            w;
            with_respect_to_first_variable,
        ),
        gfs,
    )
gradient_greens(
    gfs::_GreensFunctionVector,
    e1,
    e2,
    w;
    with_respect_to_first_variable=false,
) = sum(
    gradient_greens(gf, e1, e2, w; with_respect_to_first_variable)
    for gf in gfs
)


"""
    integral(::GreensFunction, element_1, element_2, wavenumber)

Calculates the integral of Green's function over a panel.

# Arguments
- `element_1`: first point, as a `NamedTuple` with the field `center`.
- `element_2`: integration panel, as a `NamedTuple` with the fields `center`, `normal`, `vertices`, `radius`, `normal` and `area`.
- `wavenumber`: Defines the wavenumber for waves ~~~~~. Optional for some Greens functions such as `Rankine()`.
"""
function integral end

integral(gfs::Tuple, e1, e2, w) =
    _tuple_sum(gf -> integral(gf, e1, e2, w), gfs)
integral(gfs::_GreensFunctionVector, e1, e2, w) =
    sum(integral(gf, e1, e2, w) for gf in gfs)


"""
    integral_gradient(::GreensFunction, element_1, element_2, wavenumber; with_respect_to_first_variable=false)

Calculates the integral of the gradient of Green's function over a panel.

# Arguments
- `element_1`: first point, as a `NamedTuple` with the field `center`.
- `element_2`: integration panel, as a `NamedTuple` with the fields `center`, `normal`, `vertices`, `radius`, `normal` and `area`.
- `wavenumber`: Defines the wavenumber for waves ~~~~~. Optional for some Greens functions such as `Rankine()`.
- `with_respect_to_first_variable`: A boolean flag (default is `false`). If `true`, computes the gradient with respect to `element_1`; otherwise, computes the gradient with respect to `element_2`.
"""
function integral_gradient end

integral_gradient(gfs::Tuple, e1, e2, w; with_respect_to_first_variable=false) =
    _tuple_sum(
        gf -> integral_gradient(
            gf,
            e1,
            e2,
            w;
            with_respect_to_first_variable,
        ),
        gfs,
    )
integral_gradient(
    gfs::_GreensFunctionVector,
    e1,
    e2,
    w;
    with_respect_to_first_variable=false,
) = sum(
    integral_gradient(gf, e1, e2, w; with_respect_to_first_variable)
    for gf in gfs
)


"""
    both_integral_and_integral_gradient(::GreensFunction, element_1, element_2, wavenumber; with_respect_to_first_variable=false)
Calculates both the integral of Green's function and the integral of the gradient of Green's function over a panel.
# Arguments
- `element_1`: first point, as a `NamedTuple` with the field `center`.
- `element_2`: integration panel, as a `NamedTuple` with the fields `center`, `normal`, `vertices`, `radius`, `normal` and `area`.
- `wavenumber`: Defines the wavenumber for waves ~~~~~. Optional for some Greens functions such as `Rankine()`.
- `with_respect_to_first_variable`: A boolean flag (default is `false`). If `true`, computes the gradient with respect to `element_1`; otherwise, computes the gradient with respect to `element_2`.
"""
function both_integral_and_integral_gradient(gf, element_1, element_2, wavenumber; with_respect_to_first_variable=false)
    # Default implementation, meant to be overloaded for better performance
    return (integral(gf, element_1, element_2, wavenumber),
            integral_gradient(gf, element_1, element_2, wavenumber; with_respect_to_first_variable=with_respect_to_first_variable))
end

both_integral_and_integral_gradient(
    gfs::Tuple,
    e1,
    e2,
    w;
    with_respect_to_first_variable=false,
) = _tuple_broadcast_sum(
    gf -> both_integral_and_integral_gradient(
        gf,
        e1,
        e2,
        w;
        with_respect_to_first_variable,
    ),
    gfs,
)
both_integral_and_integral_gradient(
    gfs::_GreensFunctionVector,
    e1,
    e2,
    w;
    with_respect_to_first_variable=false,
) = reduce(
    .+,
    both_integral_and_integral_gradient(
        gf,
        e1,
        e2,
        w;
        with_respect_to_first_variable,
    ) for gf in gfs
)


# Below some common tools used for several wave terms

function with_reduced_coordinates(𝒢, element_1, element_2, wavenumber)
    x  = center(element_1)
    xi = center(element_2)
    
    dis_squared_raw = (x[1] - xi[1])^2 + (x[2] - xi[2])^2
    purt = 1e-18
    if dis_squared_raw==0
        # must perturb before sqrt since sqrt(0) returns NaN dual number derivative
        dis = sqrt(dis_squared_raw + purt^2) - purt
    else
        dis = sqrt(dis_squared_raw)
    end
    r = wavenumber * dis

    # r  = wavenumber * hypot(x[1] - xi[1], x[2] - xi[2])
    # if r == 0.0
    #     r = r + 1e-18
    # end

    z  = wavenumber * (x[3] + xi[3])
    return wavenumber * 𝒢(r, z)
end

function with_reduced_coordinates_derivative(d𝒢, element_1, element_2, wavenumber; with_respect_to_first_variable=false)
    x  = center(element_1)
    xi = center(element_2)

    dis_squared_raw = (x[1] - xi[1])^2 + (x[2] - xi[2])^2
    purt = 1e-18
    if dis_squared_raw==0
        # must perturb before sqrt since sqrt(0) returns NaN dual number derivative
        dis = sqrt(dis_squared_raw + purt^2) - purt
    else
        dis = sqrt(dis_squared_raw)
    end

    r = wavenumber * dis

    # r  = wavenumber * hypot(x[1] - xi[1], x[2] - xi[2])
    # if r == 0.0
    #     r = r + 1e-18
    # end

    z  = wavenumber * (x[3] + xi[3])
    GF, dGF_dr, dGF_dz = d𝒢(r, z)
    if with_respect_to_first_variable
        if abs(r) > 1e-6
            dr_dx1 = wavenumber^2 / r * (x[1] - xi[1])
            dr_dx2 = wavenumber^2 / r * (x[2] - xi[2])
        else
            dr_dx1 = zero(x[1])
            dr_dx2 = zero(x[2])
        end
        dz_dx = wavenumber
        return wavenumber * GF, wavenumber * (zero(x) .+ (dr_dx1 * dGF_dr, dr_dx2 * dGF_dr, dz_dx * dGF_dz))
        # The zero(x) is a workaround to convert the following tuple to the same type as `x` (either Vector or SVector).
    else
        if abs(r) > 1e-6
            dr_dxi1 = wavenumber^2 / r * (xi[1] - x[1])
            dr_dxi2 = wavenumber^2 / r * (xi[2] - x[2])
        else
            dr_dxi1 = zero(x[1])
            dr_dxi2 = zero(x[2])
        end
        dz_dxi = wavenumber
        return wavenumber * GF, wavenumber * (zero(x) .+ (dr_dxi1 * dGF_dr, dr_dxi2 * dGF_dr, dz_dxi * dGF_dz))
        # The zero(x) is a workaround to convert the following tuple to the same type as `x` (either Vector or SVector).
    end
end
