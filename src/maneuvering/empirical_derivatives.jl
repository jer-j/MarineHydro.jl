raw"""
    HullParticulars(; length_pp, beam, draft, block_coefficient)

The principal dimensions every empirical maneuvering regression is built on.

`length_pp` is between perpendiculars, `beam` moulded, `draft` at even keel, and
`block_coefficient` the usual ``\nabla/(L B T)``.
"""
struct HullParticulars{T}
    length_pp::T
    beam::T
    draft::T
    block_coefficient::T
end

function HullParticulars(; length_pp::Real, beam::Real, draft::Real,
        block_coefficient::Real)
    all(>(0), (length_pp, beam, draft, block_coefficient)) ||
        throw(ArgumentError("all particulars must be positive"))
    return HullParticulars(promote(length_pp, beam, draft, block_coefficient)...)
end

raw"""
    empirical_sway_derivatives(hull; method)

Linear sway velocity derivatives ``Y'_v`` and ``N'_v`` from a named empirical or
theoretical method, in **MMG normalisation** — forces on ``\tfrac12\rho L T U^2``
and moments on ``\tfrac12\rho L^2 T U^2``.

Sign convention here follows the maneuvering literature these formulas are
quoted in, where ``Y'_v`` is reported **positive**; the hydrodynamic derivative
itself is negative. Use [`mmg_to_wang_velocity_derivatives`](@ref) to convert to
the ``\tfrac12\rho L^2`` normalisation the rest of this package uses.

Available `method` values, with the form each takes (all as collected by
Chame et al. 2025, Table 1, who in turn cite the originals):

| method | ``Y'_v`` | ``N'_v`` |
|---|---|---|
| `:low_aspect_ratio` | ``\pi T/L`` | ``\pi T/2L`` |
| `:jacobs` | ``\pi T/L + C_o`` | ``(\pi T/2L)(2x_p/L)`` |
| `:smitt` | ``1.59\,\pi T/L`` | ``1.24\,\pi T/2L`` |
| `:norrbin` | ``(\pi T/L)(1.69+0.08\,C_B B/\pi T)`` | ``(\pi T/L)(0.64-0.04\,C_B B/\pi T)`` |
| `:inoue` | ``(\pi T/L)(1+1.4\,C_B B/\pi T)`` | ``(\pi T/L)(\pi/2)`` |
| `:clarke` | ``(\pi T/L)(1+0.4\,C_B B/T)`` | ``(\pi T/L)(0.5+2.4\,T/L)`` |
| `:ho_young_lee` | ``0.4545-0.065\,C_B B/T`` | ``0.23-0.0059\,T/L`` |

`:low_aspect_ratio` is the slender-body / low-aspect-ratio wing result and is the
theoretical basis the regressions were built around; the rest add empirical
corrections fitted to model-test databases, most of them assembled decades ago.
That vintage is the point of Chame et al.: on modern hull forms these
regressions were found to over-predict ``Y'_v`` substantially.
"""
function empirical_sway_derivatives(hull::HullParticulars; method::Symbol,
        jacobs_offset::Real = 0.0, jacobs_centre::Real = 0.0)
    length_pp = hull.length_pp
    beam = hull.beam
    draft = hull.draft
    block = hull.block_coefficient
    ratio = π * draft / length_pp
    slenderness = draft / length_pp

    if method === :low_aspect_ratio
        return (Y_v = ratio, N_v = ratio / 2)
    elseif method === :jacobs
        return (Y_v = ratio + jacobs_offset,
            N_v = ratio / 2 * (2 * jacobs_centre / length_pp))
    elseif method === :smitt
        return (Y_v = 1.59 * ratio, N_v = 1.24 * ratio / 2)
    elseif method === :norrbin
        shape = block * beam / (π * draft)
        return (Y_v = ratio * (1.69 + 0.08 * shape),
            N_v = ratio * (0.64 - 0.04 * shape))
    elseif method === :inoue
        return (Y_v = ratio * (1 + 1.4 * block * beam / (π * draft)),
            N_v = ratio * (π / 2))
    elseif method === :clarke
        return (Y_v = ratio * (1 + 0.4 * block * beam / draft),
            N_v = ratio * (0.5 + 2.4 * slenderness))
    elseif method === :ho_young_lee
        return (Y_v = 0.4545 - 0.065 * block * beam / draft,
            N_v = 0.23 - 0.0059 * slenderness)
    end
    throw(ArgumentError("unknown method $(method)"))
end

"""
    empirical_methods()

Names accepted by [`empirical_sway_derivatives`](@ref), in the order they are
usually tabulated.
"""
empirical_methods() = (:low_aspect_ratio, :jacobs, :smitt, :norrbin, :inoue,
    :clarke, :ho_young_lee)

raw"""
    clarke_full_derivatives(hull)

The complete Clarke, Gedling and Hine (1982) set — both acceleration and
velocity derivatives — in MMG normalisation.

Every derivative takes the same form,

```math
(\,\cdot\,)'=-\pi\left(\frac{T}{L}\right)^{2}
\left[a_0+a_1\frac{B}{L}+a_2\frac{B}{T}+a_3\frac{T}{L}
+a_4 C_B\frac{B}{T}+a_5\left(\frac{B}{L}\right)^{2}
+a_6\left(\frac{B}{T}\right)^{2}\right],
```

with the coefficients of Gokarn, *A Study of Ship Manoeuvrability*, Table 11.2.
Returned with the sign the equation gives, so ``Y_v`` and ``N_r`` come out
negative.
"""
function clarke_full_derivatives(hull::HullParticulars)
    beam_length = hull.beam / hull.length_pp
    beam_draft = hull.beam / hull.draft
    draft_length = hull.draft / hull.length_pp
    block = hull.block_coefficient
    prefactor = -π * draft_length^2

    coefficients = (Y_vdot = (1.0, 0.0, 0.0, 0.0, 0.16, -5.1, 0.0),
        Y_v = (1.0, 0.0, 0.0, 0.0, 0.4, 0.0, 0.0),
        N_vdot = (0.0, 1.1, -0.041, 0.0, 0.0, 0.0, 0.0),
        N_v = (0.5, 0.0, 0.0, 2.4, 0.0, 0.0, 0.0),
        Y_rdot = (0.0, 0.67, 0.0, 0.0, 0.0, 0.0, -0.0033),
        Y_r = (-0.5, 0.0, -0.08, 0.0, 0.0, 2.2, 0.0),
        N_rdot = (0.0833, -0.33, 0.0, 0.0, 0.017, 0.0, 0.0),
        N_r = (0.25, -0.56, 0.0, 0.0, 0.039, 0.0, 0.0))

    basis = (1.0, beam_length, beam_draft, draft_length, block * beam_draft,
        beam_length^2, beam_draft^2)
    return map(row -> prefactor * sum(row .* basis), coefficients)
end

raw"""
    hirano_takashina_derivatives(hull)

The low-aspect-ratio set of Hirano and Takashina (2010), as given by Gokarn,
in MMG normalisation with ``\Lambda=2T/L`` the effective aspect ratio of the
hull and its free-surface image.

```math
Y'_v=-\tfrac12\pi\Lambda+1.4\,C_B\frac{B}{L},\quad
N'_v=-\Lambda,\quad
Y'_r=\tfrac14\pi\Lambda,\quad
N'_r=-\left(0.54\Lambda-\Lambda^{2}\right)
```
"""
function hirano_takashina_derivatives(hull::HullParticulars)
    aspect = 2 * hull.draft / hull.length_pp
    return (Y_v = -0.5π * aspect +
                  1.4 * hull.block_coefficient * hull.beam / hull.length_pp,
        N_v = -aspect, Y_r = 0.25π * aspect,
        N_r = -(0.54 * aspect - aspect^2))
end
