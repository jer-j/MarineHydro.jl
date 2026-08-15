# KVLCC2 maneuvering references

## Source and role

`reference.toml` transcribes Tables 4, 7, and 8 from:

In-Tae Kim et al., "Estimation of the manoeuvrability of the KVLCC2 in calm
water using free running simulation based on CFD," *International Journal of
Naval Architecture and Ocean Engineering* 13 (2021), 466-477,
[doi:10.1016/j.ijnaoe.2021.05.004](https://doi.org/10.1016/j.ijnaoe.2021.05.004).

The paper is not a captive-test derivative report. It combines free-running
RANS simulations, SIMMAN 2008 experimental maneuver results, and a
system-based MMG calculation. The different data are kept separate here:

- Table 4 supplies model-test-derived MMG hull coefficients for an immediate
  secondary comparison with the bare-hull potential-flow calculation.
- Tables 7 and 8 supply experimental zig-zag overshoot angles and turning
  indices. They are future end-to-end targets for a simulator that includes
  hull, propeller, rudder, and their interactions.

## Linear coefficient conversion

The standard MMG definitions use

```math
Y_H=\frac12\rho L_{pp}TU^2Y_H',\qquad
N_H=\frac12\rho L_{pp}^2TU^2N_H',
```

with $v'=v/U$ and $r'=rL_{pp}/U$. Wang uses one additional power of
$L_{pp}$ instead of $T$ in each corresponding derivative scale. Therefore

```math
\left(Y_v',Y_r',N_v',N_r'\right)_W
=\frac{T}{L_{pp}}
\left(Y_v',Y_r',N_v',N_r'\right)_{MMG}.
```

For KVLCC2, $T/L_{pp}=20.8/320=0.065$, giving

| Coefficient | Table 4 MMG | Wang normalization |
|---|---:|---:|
| $Y_v'$ | -0.315 | -0.020475 |
| $Y_r'$ | 0.083 | 0.005395 |
| $N_v'$ | -0.137 | -0.008905 |
| $N_r'$ | -0.049 | -0.003185 |

The paper uses $y$ to starboard and $z$ downward. MarineHydro uses $y$
to port and $z$ upward. This is a 180-degree rotation about $x$: each
generalized velocity and its conjugate load change sign together, leaving the
four first-order derivative signs unchanged.

## Scientific interpretation

The Table 4 coefficients contain viscous and empirical MMG content. They are
not expected to equal whole-hull, inviscid double-body coefficients. The
relevant Wang comparison uses Schmitz truncation at the aft end of the
maximum-area region. `gothenburg_maximum_section` derives that cutoff from the
official structured KVLCC2 surface rather than selecting a station by hand.

Run the comparison from
[`notebooks/02_kvlcc2_viscous_correction.ipynb`](../../notebooks/02_kvlcc2_viscous_correction.ipynb),
which reports whole-hull, Schmitz-truncated and viscously corrected velocity
derivatives against the values below over a four-mesh sequence, and writes its
figures to `validation/gothenburg2010/results/`.

## Current comparison

On the 480-panel mesh, the automatic criterion selects
$x_{cut}/L_{pp}=-0.095008$ at the aft end of the near-maximum-area region.
The discrete maximum itself is at $x/L_{pp}=-0.077759$. The cutoff retains
256 of the 480 panels in the velocity-dependent surface integrals.

The independent unit-sway and unit-yaw calculation, with the current
`base_flow = :double_body` and `acceleration_formulation = :direct` defaults,
gives:

| Coefficient | Schmitz BEM | Table 4 converted | Relative difference |
|---|---:|---:|---:|
| $Y_v'$ | -0.018794 | -0.020475 | 8.21% |
| $Y_r'$ | -0.001299 | +0.005395 | 124.08%, opposite sign |
| $N_v'$ | -0.007785 | -0.008905 | 12.58% |
| $N_r'$ | -0.002078 | -0.003185 | 34.76% |

Clarke's regressions, which is what Wang et al. actually recommend for the
rotational pair, give $Y_r'=+0.003817$ and $N_r'=-0.002549$ from the reference
$Y_v'$ — correct signs and within 30%. See `clarke_rotational_derivatives`.

The truncation is essential. Whole-hull integration gives $Y_v'=+0.001111$ and
$N_v'=-0.014839$; Schmitz integration changes these to $-0.018794$ and
$-0.007785$, much closer to the model-test-derived MMG coefficients. This
supports the implementation of Wang's stated integration domain, but it does
not validate $Y_r'$ and it is not a substitute for raw PMM data.

These numbers moved when the incomplete steady linearization was corrected (see
`docs/WANG_MANEUVERING_PLAN.md`). Under Wang's original uniform-stream form the
same mesh gave $Y_v'=-0.021803$ and $N_v'=-0.009082$, i.e. 6.48% and 1.98% from
the references — *better* agreement than the more correct linearization
produces. That earlier agreement was partly a cancellation between the
linearization error and the error in the semi-empirical stern cut, which is a
reason to treat the truncated velocity derivatives as calibrated rather than
predicted.
