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

The code implements Stages 1 and 2 through:

- `solve_wang_maneuvering`
- `evaluate_indirect_potential_gradient`
- `wang_hydrodynamic_derivatives`
- `wang_stern_mask`
- `nondimensionalize_wang_derivatives`
- `wang_restricted_water_elevation` for Wang et al.'s sectionwise continuity
  and Bernoulli equations

The current implementation also includes:

- `solve_potential_flow_maneuvering` for independent sway and yaw potentials
- `read_tecplot_structured_surfaces` and `read_gothenburg2010_mesh`
- `structured_sectional_area_curve` and `gothenburg_maximum_section` for an
  automatic Schmitz cutoff derived from the official structured hull surface
- `MMGLinearHullDerivatives` and `mmg_to_wang_velocity_derivatives` for an
  explicit standard-MMG to Wang normalization conversion
- checksum-verified KVLCC2 and KCS geometry acquisition
- KVLCC2 geometry, derivative-convergence, and MMG-comparison drivers
- 1,984-panel surge-flow orientation and surface-streamline output
- a cited KVLCC2 reference dataset containing the linear hull coefficients
  and experimental maneuver indices transcribed from Kim et al. (2021)

The tests exercise the low-level quadrature, potential-gradient trace, both
boundary residuals, acceleration reciprocity, speed scaling, truncation mask,
nondimensionalization, coordinate transformation, mirroring, and panel
orientation. The mesh storage remains generic over concrete array and scalar
types, including `ForwardDiff.Dual`, while retaining enough type information
for specialization. JET optimization analysis reports zero findings for the
boundary conditions, derivative integrals, benchmark importer, sectional-area
cutoff, normalization conversion, complete BEM solves, restricted-water
elevation, stagnation search, and surface-streamline tracer. The package-wide
test suite passes 1,804 tests, including the existing Zygote and ForwardDiff
coverage.

The official-geometry KVLCC2 convergence outputs are bit-for-bit unchanged by
the inference improvements. On the 480-panel mesh, automatic Schmitz
truncation gives $Y_v'=-0.021803$ and $N_v'=-0.009082$, within 6.48% and
1.98% of the converted model-test-derived MMG values. The independent-yaw
$Y_r'$ has the opposite sign, so that coefficient is not validated. Raw
SIMMAN captive-test ingestion, a finer mesh study, and authoritative
Mariner/Tokyo Maru geometry remain the next ship-data milestones.
