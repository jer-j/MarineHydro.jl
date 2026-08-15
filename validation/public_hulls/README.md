# Public-hull derivative and added-mass study

This validation stage applies MarineHydro's unrestricted, infinite-depth
double-body boundary-element method to five publicly reproducible hulls:

- KCS and KVLCC2 from the checksum-verified NMRI Gothenburg workshop surfaces
- DTMB 5415 from NMRI's formatted Plot3D `5415_static.net` surface
- DTC from the pinned OpenFOAM `DTC-scaled.stl.gz` resource
- the analytic Wigley hull used in Ayub's 2021 Kyushu University thesis

The calculation solves independent unit-sway and unit-yaw potentials. Wang's
strict approximation, $\phi_r\approx x\phi_v$, is evaluated from the same sway
solution. Acceleration derivatives use every wetted panel. The configured
Schmitz cutoff is applied only to the convective velocity-derivative integrals.

## Reproduce

Fetch all external geometry with pinned SHA-256 checksums:

```sh
validation/public_hulls/fetch_geometry.sh
```

Run the standard comparison, which remains below 500 panels per hull:

```sh
julia --project=. validation/public_hulls/public_hull_validation.jl
```

Run the higher-resolution study with:

```sh
julia --project=. \
  validation/public_hulls/public_hull_validation.jl --fine
```

The Julia/IJulia notebook
[`notebooks/public_hull_derivatives.ipynb`](../../notebooks/public_hull_derivatives.ipynb)
executes the same functions and produces a separate three-dimensional view of
every hull, static-drift comparisons, and added-mass comparisons. Generated CSV
and PNG files are ignored under `results/`. No free-surface elevation or plane
is rendered.

## Definitions

The added-mass matrix follows

```math
A_{22}=-Y_{\dot v},\qquad
A_{26}=-Y_{\dot r},\qquad
A_{62}=-N_{\dot v},\qquad
A_{66}=-N_{\dot r}.
```

The comparison reports both displacement normalization and standard MMG
normalization:

```math
C_{A22}=\frac{A_{22}}{\rho\nabla},\qquad
C_{A66}=\frac{A_{66}}{\rho\nabla L_{pp}^2},
```

```math
m_y'=\frac{A_{22}}{\tfrac12\rho L_{pp}^2T},\qquad
J_z'=\frac{A_{66}}{\tfrac12\rho L_{pp}^4T}.
```

With $\beta\simeq-v/U$ in MarineHydro's $y$-port convention, the reported
static-drift slopes are

```math
Y_\beta'=-\frac{L_{pp}}{T}Y_{v,W}',\qquad
N_\beta'=-\frac{L_{pp}}{T}N_{v,W}'.
```

## Standard-resolution checkpoint

The implementation checkpoint below uses 288 to 479 panels. These values are
intended to catch regressions and establish trends, not to replace the `--fine`
convergence study. The complete numeric checkpoint is tracked in
[`standard_checkpoint.csv`](standard_checkpoint.csv).

| Hull | Panels | $m_y'$ | Reference | $J_z'$ | Reference | $Y_\beta'$ | Reference | $N_\beta'$ | Reference |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| KCS | 364 | 0.1535 | 0.1552 | 0.00879 | 0.00927 | 0.3355 | 0.2252 | 0.1228 | 0.1111 |
| KVLCC2 | 364 | 0.2395 | 0.2230 | 0.01267 | 0.0110 | 0.3358 | 0.3150 | 0.1407 | 0.1370 |
| DTMB 5415 | 462 | 0.1145 | unavailable | 0.00767 | unavailable | 0.1159 | 0.1160 | 0.0319 | 0.0640 |
| DTC | 479 | 0.1274 | unavailable | 0.00695 | unavailable | 0.2255 | unavailable | 0.0775 | unavailable |
| Wigley | 288 | 0.1670 | unavailable | 0.01282 | unavailable | 0.4710 | 0.2583 | 0.2032 | 0.1248 |

Additional displacement-normalized sway results are
$A_{22}/(\rho\nabla)=0.8423$, $0.8157$, $0.8387$, $0.6712$, and $1.8788$ for
KCS, KVLCC2, DTMB 5415, DTC, and Wigley, respectively. The approximate
published DTC deep-water value is $0.70$. At this resolution, volume errors
are $-1.02\%$, $-1.02\%$, $-1.05\%$, $-3.52\%$, and $-0.81\%$. Linear-system
residuals are below $2.3\times10^{-15}$, and acceleration reciprocity errors
remain below $0.50\%$.

The KCS and KVLCC2 added masses are encouraging. The Schmitz-truncated KVLCC2
drift slopes are also close to the model-test-derived MMG values. The DTMB yaw
slope and both Wigley empirical comparisons show that inviscid pressure plus a
single stern cutoff is not a complete maneuvering model. Viscous cross-flow,
finite-Froude-number effects, cutoff sensitivity, and raw captive-test
uncertainty must remain separate validation layers.

## Public sources and reference status

- NMRI DTMB 5415 geometry and condition page:
  <https://www.nmri.go.jp/archives/institutes/fluid_performance_evaluation/cfd_rd/cfdws05/gothenburg2000/5415/5415_g%26c.htm>
- OpenFOAM DTC geometry pinned at commit `e11dbc66`:
  <https://github.com/OpenFOAM/OpenFOAM-dev/blob/e11dbc66d29b34031e9cc6335ac20ae2d1487af2/tutorials/resources/geometry/DTC-scaled.stl.gz>
- Ayub, *Influence of Hydrodynamic Derivatives on Ship Manoeuvring Prediction
  and Application of SQCM to Ship Hull Forces*, Kyushu University, 2021:
  <https://hdl.handle.net/2324/4496053>
- Kim et al., KVLCC2 MMG coefficients and maneuvering experiments, 2021:
  <https://doi.org/10.1016/j.ijnaoe.2021.05.004>
- DTMB 5415 EFD derivative comparison, 2026:
  <https://doi.org/10.3390/jmse14141274>
- DTC deep-water added-mass context, 2025:
  <https://doi.org/10.1080/09377255.2025.2513760>

The KCS and KVLCC2 added-mass comparators are MMG values, not direct
zero-frequency oscillation measurements. The DTC $A_{22}/(\rho\nabla)$ value is
reported approximately in the cited paper. DTMB values are EFD static-drift
slopes at $F_n=0.41$, not added masses. Wigley values are Kijima empirical
estimates reproduced by Ayub. The CSV retains `NaN` when no authoritative
numeric comparator has been transcribed.
