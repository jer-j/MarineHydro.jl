# Gothenburg and SIMMAN validation data

## Scope

Gothenburg 2010 used the original KVLCC2, KCS, and DTMB 5415 hulls. Its 18
cases cover resistance, sinkage and trim, free-surface elevation, local flow,
self-propulsion, seakeeping, and roll decay. It does not publish captive-test
sway-yaw hydrodynamic derivatives.

This leads to two distinct validation tracks:

1. Use the Gothenburg hull definitions for mesh parsing, coordinate signs,
   surface area, displacement, and spatial convergence.
2. Use SIMMAN captive-motion measurements on the same KVLCC2 hull for
   maneuvering-derivative validation. The public SIMMAN pages give geometry,
   conditions, and test plans. The numerical model-test files require access
   requested from the SIMMAN maintainers and are therefore not redistributed
   here.

The relevant primary sources are:

- [Gothenburg 2010 workshop record](https://research.chalmers.se/en/publication/131971)
- [Gothenburg 2010 proceedings](https://doi.org/10.1007/978-94-007-7189-5)
- [NMRI KVLCC2 geometry and conditions](https://www.nmri.go.jp/archives/institutes/fluid_performance_evaluation/cfd_rd/cfdws05/gothenburg2000/KVLCC/kvlcc_g%26c.htm)
- [SIMMAN KVLCC2 geometry and conditions](https://simman2014.dk/ship-data/moeri-kvlcc2-tanker/geometry-and-conditions-moeri-kvlcc2-tanker/)
- [SIMMAN KVLCC2 test overview](https://simman2014.dk/ship-data/moeri-kvlcc2-tanker/overview-of-tests-moeri-kvlcc2-tanker/)
- [SIMMAN deep-water bare-hull PMM test specification](https://simman2014.dk/wp-content/uploads/2015/07/KVLCC2-PMM-deep-BARE-HULL.pdf)
- [Yasukawa and Yoshimura MMG standard method](https://doi.org/10.1007/s00773-014-0293-y)

Gothenburg 2010 returned to the original KVLCC2 definition after the faired
KVLCC2M variant had been used at Tokyo 2005. The NMRI archive supplies that
original KVLCC2 geometry and the same KCS definition, so the archived
structured surfaces are suitable geometric inputs for the 2010 hulls.

## Reference values

| Hull | $L_{pp}$ (m) | $B$ (m) | $T$ (m) | $\nabla$ (m³) | $S$ (m²) | $C_B$ | $F_n$ |
|---|---:|---:|---:|---:|---:|---:|---:|
| KVLCC2 | 320.0 | 58.0 | 20.8 | 312622 | 27194 | 0.8098 | 0.142 |
| KCS | 230.0 | 32.2 | 10.8 | 52030 | 9424 | 0.6505 | 0.260 |
| DTMB 5415 | 142.0 | 19.06 | 6.15 | 8424.4 | 2972.6 | 0.507 | 0.248 or 0.413 |

The KVLCC2 surface-area value excludes the rudder. SIMMAN reports
$S=27467\ \mathrm{m^2}$ when the rudder is included.

## Reproducible geometry acquisition

From the repository root, run:

```sh
validation/gothenburg2010/fetch_geometry.sh
```

The script downloads the official KVLCC2 and KCS structured-surface archives,
verifies fixed SHA-256 checksums, and extracts them into the ignored `data/`
directory. The source files are nondimensionalized by $L_{pp}$ and use
$x$ positive aft and $y$ positive starboard. `read_gothenburg2010_mesh`
converts them to MarineHydro's $x$-forward, $y$-port convention:

```math
(x,y,z)_{MH}=L_{pp}(-x,-y,z)_{G}.
```

## Running the checks

Everything that used to live in `kvlcc2_convergence.jl`,
`kvlcc2_boundary_layer_validation.jl` and `kvlcc2_flow_visualization.jl` is now
in one executable notebook:

[`notebooks/02_kvlcc2_viscous_correction.ipynb`](../../notebooks/02_kvlcc2_viscous_correction.ipynb)

It covers, in order:

- **Geometry convergence** — wetted surface and displaced volume against the
  published values, from the coarsest resampling up to the native structured
  grid.
- **Double-body surge flow** — a 1,984-panel solve rendered as full-hull, bow
  and stern three-dimensional views, marking the minimum-speed bow panel and
  tracing surface streamlines. In MarineHydro coordinates the bow is positive
  $x$ and the incoming water velocity is toward negative $x$. No free-surface
  elevation is drawn: this is the zero-frequency double-body limit, and a wave
  pattern would require a steady finite-Froude-number free-surface solve that
  the model does not contain.
- **Boundary-layer verification** — the Head closure against the turbulent
  flat-plate 1/5-power laws, and the quasi-3D strip march against the ITTC-1957
  friction line at $Re = 4.6\times10^{6}$. The predicted $C_F$ is $0.947$–$0.950$
  times ITTC-57 and is mesh-independent across a sevenfold change in panel
  count.
- **Derivative comparison** — whole-hull inviscid, Schmitz-truncated, Head shear,
  Head displacement-pressure and combined results reported as separate rows, so
  the Schmitz and physical boundary-layer corrections are never silently counted
  twice. The reference is the KVLCC2 MMG coefficient set transcribed from Kim et
  al. (2021); see the
  [KVLCC2 maneuvering reference](../kvlcc2_maneuvering/README.md) for the
  normalization, coordinate mapping and limitations.

Figures and CSVs are written to the ignored `results/` directory. The
experimental local-flow conditions and measurement locations from Case 1.1a are
transcribed in `local_flow_reference.toml`. See
[`docs/INTEGRAL_BOUNDARY_LAYER.md`](../../docs/INTEGRAL_BOUNDARY_LAYER.md) for
the equations, assumptions and validation boundary.

The notebook uses the reproducible one-way correction. A relaxed
FlightStream-style viscous–inviscid iteration is available by passing
`coupling_iterations=8` to `viscous_maneuvering_correction`; the result records
the selected linearization, the coupling residual, the BEM residuals and the
integrated differentiated transpiration flux.

On the complete public KVLCC2 structured grid, the current importer gives
$\nabla/L_{pp}^3=0.00954121$, an error of $+0.00775\%$ relative to the
published displacement. The faceted surface has $S/L_{pp}^2=0.270510$, which
is $1.86\%$ above the published smooth-hull area. KCS exhibits the same
one-sided faceting trend, while its displacement error is only $+0.0294\%$.
Surface-area and volume convergence are therefore reported independently.

## Derivative-data status

The highest-value comparison is the deep-water, bare-hull INSEAN PMM test for
KVLCC2 at scale 1:100, $L_{pp}=3.2\ \mathrm{m}$, and $F_n=0.142$. The SIMMAN
site requires a named access request for the raw FTP model data. This
repository will accept those files through a separate parser once legitimate
access is available. Published MMG coefficient sets are useful secondary
comparators but are not substitutes for the raw captive-test time series.

### Secondary added-mass check

Yasukawa and Yoshimura give $m_y'=0.223$ and $J_z'=0.011$ for KVLCC2, using
$\tfrac12\rho L_{pp}^2T$ and $\tfrac12\rho L_{pp}^4T$ as the mass and yaw
inertia scales. Their force convention contains $-m_y\dot v$ and
$-J_z\dot r$. Conversion to Wang's scales is therefore

```math
Y_{\dot v,W}'=-\frac{T}{L_{pp}}m_y'=-0.014495,
\qquad
N_{\dot r,W}'=-\frac{T}{L_{pp}}J_z'=-0.000715.
```

On the 480-panel mesh, `solve_potential_flow_maneuvering` gives
$Y_{\dot v}'=-0.015104$ and $N_{\dot r}'=-0.000808$, differences of $4.20\%$
and $12.97\%$ from those MMG values. The strict Wang yaw approximation gives
$N_{\dot r}'=-0.000972$, a $35.9\%$ difference. This is encouraging but is not
an experimental validation: the MMG paper states that its added-mass values
were estimated from Motora's empirical charts, and the BEM sequence is still
coarse.
