# Wang maneuvering derivatives for MarineHydro.jl

## Objective

Add a validated implementation of the low-frequency method of Wang,
Bakountouzis, and Katory for predicting linear sway and yaw hydrodynamic
derivatives from a boundary-element velocity potential.

The first milestone is restricted to unrestricted horizontal boundaries and
infinite water depth. It uses the zero-frequency rigid-free-surface, or
double-body, approximation. Forward speed enters the linearized Bernoulli
equation but not the boundary-integral operator.

## Governing problem

The unit-sway potential satisfies

```math
\nabla^2\phi_2=0,
\qquad
\frac{\partial\phi_2}{\partial n}=n_y\quad\text{on }S_h,
\qquad
\frac{\partial\phi_2}{\partial z}=0\quad\text{on }z=0.
```

The default Green function is

```math
G_0(\boldsymbol{x},\boldsymbol{\xi})
=
\frac{1}{|\boldsymbol{x}-\boldsymbol{\xi}|}
+
\frac{1}{|\boldsymbol{x}-\boldsymbol{\xi}^{*}|},
```

in MarineHydro's internal normalization. The reflected source has the same
sign, so the free-surface Neumann condition is satisfied exactly.

Wang et al. also calculate a mean water-level change in restricted channels
from sectionwise continuity and Bernoulli equations. MarineHydro exposes this
as `wang_restricted_water_elevation`. In the present unrestricted limit,

```math
A_0,W\rightarrow\infty,
\qquad U_1(x)\rightarrow U_0,
\qquad \zeta(x)\rightarrow0.
```

Consequently, the zero-frequency double-body solution has a flat rigid free
surface. A nonzero Kelvin-wave elevation would require the deferred steady
finite-Froude-number free-surface problem and cannot be reconstructed from the
current reflected-Rankine solution.

The strict Wang derivatives are

```math
\begin{aligned}
Y_{\dot v} &= \rho\int_{S_h}\phi_2 n_y\,\mathrm{d}S,
&Y_v &= -\rho U\int_{S_h^*}\phi_{2,x} n_y\,\mathrm{d}S,\\
Y_{\dot r} &= \rho\int_{S_h}x\phi_2 n_y\,\mathrm{d}S,
&Y_r &= -\rho U\int_{S_h^*}(\phi_2+x\phi_{2,x})n_y\,\mathrm{d}S,\\
N_{\dot v} &= \rho\int_{S_h}x\phi_2 n_y\,\mathrm{d}S,
&N_v &= -\rho U\int_{S_h^*}x\phi_{2,x}n_y\,\mathrm{d}S,\\
N_{\dot r} &= \rho\int_{S_h}x^2\phi_2 n_y\,\mathrm{d}S,
&N_r &= -\rho U\int_{S_h^*}x(\phi_2+x\phi_{2,x})n_y\,\mathrm{d}S.
\end{aligned}
```

## Development stages

### Stage 1: zero-frequency numerical kernel

- Reuse `Rankine()` and `RankineReflected()`.
- Solve the unit-sway problem with the indirect formulation.
- Recover the fluid-side potential-gradient trace from the single-layer
  source strengths.
- Verify each gradient component against the existing directional derivative
  matrices.
- Verify that the reconstructed normal derivative recovers `n_y`.

Exit criteria:

- Relative boundary-condition residual below `1e-10` on canonical meshes.
- Gradient components agree with the existing matrix assembly to roundoff.

### Stage 2: strict Wang extraction

- Implement the eight dimensional derivatives.
- Support an explicit velocity-integration mask.
- Provide a helper for longitudinal Schmitz stern truncation.
- Implement the nondimensionalization from Wang et al.
- Keep the potential-flow whole-hull result available as the baseline.

Exit criteria:

- Exact discrete identity `Y_rdot == N_vdot`.
- Acceleration derivatives independent of forward speed.
- Velocity derivatives linear in forward speed.
- Correct dimensional scaling under uniform geometric rescaling.

### Stage 3: canonical validation

- Add mesh-refinement studies for a surface-piercing ellipsoid and a
  wall-sided canonical hull.
- Compare low-frequency acceleration derivatives with MarineHydro radiation
  added mass as frequency tends to zero.
- Compare the ellipsoid acceleration terms with analytical added-mass values.
- Check symmetry, reference-point transformations, and fore-aft cancellation.
- Quantify sensitivity to the stern truncation station.

Exit criteria:

- Demonstrated asymptotic convergence under mesh refinement.
- Cross-method agreement within the documented discretization uncertainty.
- All coordinate and force signs verified by independent identities.

### Stage 4: exact three-dimensional yaw extension

The strict Wang approximation represents yaw with `x * phi_2`. The implemented
general formulation solves a separate unit-yaw potential with

```math
\frac{\partial\phi_r}{\partial n}=x n_y-y n_x.
```

Generalized force modes form the full two-by-two sway-yaw acceleration and
convective-velocity matrix. The strict Wang result remains available for
reproduction of the paper.

### Stage 5: ship-data validation

- Reproduce an open-water Mariner or comparable benchmark.
- Document the geometric definition of the maximum-area truncation station.
- Compare whole-hull, Wang-Schmitz, and empirical Clarke values.
- Propagate mesh and truncation uncertainty into the reported derivatives.

## Deferred scope

The following capabilities are intentionally deferred:

- Finite-depth Green functions.
- Channel sidewalls and restricted-water blockage.
- Coupling the implemented sectionwise squat and water-level correction back
  into the wetted geometry and derivative integrals.
- Steady wave-making at finite Froude number.
- Nonlinear derivatives such as `Y_vv` or `Y_vvv`.
- Explicit viscous wakes and Kutta conditions.
- Fast multipole or hierarchical-matrix acceleration.

## Current implementation status

The code implements Stages 1 and 2 and the analytical-ellipsoid portion of
Stage 3 through:

- `solve_wang_maneuvering`
- `evaluate_indirect_potential_gradient`
- `wang_hydrodynamic_derivatives`
- `wang_stern_mask`
- `nondimensionalize_wang_derivatives`
- `wang_restricted_water_elevation` for Wang et al.'s sectionwise continuity
  and Bernoulli equations
- `ellipsoid_potential_coefficients` and `ellipsoid_added_mass` for the
  classical triaxial-ellipsoid reference solution
- `surface_piercing_ellipsoid_panel_grid` for a differentiable structured
  lower-half ellipsoid mesh

The current implementation also includes:

- `solve_potential_flow_maneuvering` for independent sway and yaw potentials
- `read_tecplot_structured_surfaces` and `read_gothenburg2010_mesh`
- `read_plot3d_structured_surfaces`, `structured_surface_waterline_mesh`, and
  `read_ascii_stl_waterline_mesh` for public DTMB 5415 and DTC geometry
- `wigley_hull_panel_grid` for an AD-compatible analytic Wigley benchmark
- `structured_sectional_area_curve` and `gothenburg_maximum_section` for an
  automatic Schmitz cutoff derived from the official structured hull surface
- `MMGLinearHullDerivatives` and `mmg_to_wang_velocity_derivatives` for an
  explicit standard-MMG to Wang normalization conversion
- checksum-verified KVLCC2 and KCS geometry acquisition
- KVLCC2 geometry, derivative-convergence, and MMG-comparison drivers
- 1,984-panel surge-flow orientation and surface-streamline output
- a cited KVLCC2 reference dataset containing the linear hull coefficients
  and experimental maneuver indices transcribed from Kim et al. (2021)
- an ellipsoid mesh-refinement driver and a cross-method comparison against
  the low-frequency limit of the independent radiation solver
- native Julia/IJulia notebooks for the ellipsoid convergence plots and the
  1,984-panel KVLCC2 bow, stern, and stagnation views
- a Julia/IJulia public-hull notebook covering KCS, KVLCC2, DTMB 5415, DTC,
  and Wigley geometry, derivatives, and added masses

The tests exercise the low-level quadrature, potential-gradient trace, both
boundary residuals, acceleration reciprocity, speed scaling, truncation mask,
nondimensionalization, coordinate transformation, mirroring, and panel
orientation. The public-hull additions test formatted Plot3D parsing, exact
waterline clipping, clustered ASCII STL import, Wigley volume convergence, and
differentiated Wigley geometry. The mesh storage remains generic over concrete
array and scalar types, including `ForwardDiff.Dual`, while retaining enough
type information for specialization. JET optimization analysis reports zero
findings for the maneuvering boundary conditions, derivative integrals,
Gothenburg importer, sectional-area cutoff, normalization conversion, complete
BEM solves, restricted-water correction, stagnation search, surface-streamline
tracer, analytical ellipsoid coefficients, ellipsoid added mass, and
structured ellipsoid mesh.

The official-geometry KVLCC2 convergence outputs are bit-for-bit unchanged by
the inference improvements. On the 480-panel mesh, automatic Schmitz
truncation gives $Y_v'=-0.021803$ and $N_v'=-0.009082$, within 6.48% and
1.98% of the converted model-test-derived MMG values. The independent-yaw
$Y_r'$ has the opposite sign, so that coefficient is not validated. Raw
SIMMAN captive-test ingestion, a finer mesh study, and authoritative
Mariner/Tokyo Maru geometry remain the next ship-data milestones.

For the canonical ellipsoid with ``(a,b,c)=(3,1,0.8)`` m, the analytical
lower-half sway added mass is ``3340.618`` kg at ``\rho=1000`` kg/m³. The
default constant-panel Wang sequence reduces its error from ``+7.87\%`` at 32
panels to ``+5.84\%`` at 200 panels. An independent 128-panel radiation
calculation gives ``3321.805`` kg at ``\omega=0.01`` rad/s, a ``-0.563\%``
analytical error, with damping at roundoff scale. This completes the
analytical and low-frequency portions of Stage 3. The wall-sided canonical
hull, systematic stern-cutoff sensitivity, and raw captive-test ingestion
remain open validation gates.
