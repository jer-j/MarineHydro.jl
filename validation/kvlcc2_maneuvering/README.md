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

Run the fine comparison with:

```sh
julia --project=. validation/gothenburg2010/kvlcc2_convergence.jl --fine
```

This writes the long-form
`validation/gothenburg2010/results/kvlcc2_mmg_velocity_comparison.csv` in
addition to the existing geometry and BEM convergence files.

## Current comparison

On the 480-panel mesh, the automatic criterion selects
$x_{cut}/L_{pp}=-0.095008$ at the aft end of the near-maximum-area region.
The discrete maximum itself is at $x/L_{pp}=-0.077759$. The cutoff retains
256 of the 480 panels in the velocity-dependent surface integrals.

The independent unit-sway and unit-yaw calculation gives:

| Coefficient | Schmitz BEM | Table 4 converted | Relative difference |
|---|---:|---:|---:|
| $Y_v'$ | -0.021803 | -0.020475 | 6.48% |
| $Y_r'$ | -0.002195 | +0.005395 | 140.69%, opposite sign |
| $N_v'$ | -0.009082 | -0.008905 | 1.98% |
| $N_r'$ | -0.002533 | -0.003185 | 20.48% |

For comparison, the strict Wang yaw approximation gives $Y_r'=+0.000830$.
It has the reference sign but is 84.62% smaller in magnitude. Its other three
relative differences are 6.48%, 6.70%, and 25.72%, respectively.

The truncation is essential. In the independent-potential calculation,
whole-hull integration gives $Y_v'=+0.000686$ and $N_v'=-0.017252$;
Schmitz integration changes these to $-0.021803$ and $-0.009082$. The
latter values are much closer to the model-test-derived MMG coefficients.
This result supports the implementation of Wang's stated integration domain,
but it does not validate $Y_r'$ and it is not a substitute for raw PMM data.
