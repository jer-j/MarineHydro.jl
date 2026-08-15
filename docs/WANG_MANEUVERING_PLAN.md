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

- `solve_potential_flow_maneuvering` for independent sway and yaw potentials,
  with `acceleration_formulation` and `base_flow` controlling the two
  corrections described under "Corrections applied after verification"
- `steady_pressure_derivative` for the complete body-frame linearization of the
  steady pressure with respect to sway velocity or yaw rate
- `clarke_rotational_derivatives` for Wang's equations (31) and (32), and
  `rotational_velocity_mask` for a separate rotational truncation station
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
- an ellipsoid mesh-refinement study and a cross-method comparison against
  the low-frequency limit of the independent radiation solver
- three native Julia/IJulia notebooks that replace the former standalone
  validation drivers and carry the theory alongside the calculation:
  `notebooks/01_maneuvering_theory.ipynb`,
  `notebooks/02_kvlcc2_viscous_correction.ipynb`, and
  `notebooks/03_public_hull_derivatives.ipynb`

The tests exercise the low-level quadrature, potential-gradient trace, both
boundary residuals (see the note on their limited meaning at the end of this
document), acceleration reciprocity, speed scaling, truncation mask,
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

On the 480-panel official-geometry KVLCC2 mesh, automatic Schmitz truncation
gives $Y_v'=-0.018794$ and $N_v'=-0.007785$, within 8.21% and 12.58% of the
converted model-test-derived MMG values. Directly integrated $Y_r'$ still has
the opposite sign, because it needs the separate rotational truncation station
or the Clarke regressions of §3 below rather than the sway cut. Raw SIMMAN
captive-test ingestion, a finer mesh study, and authoritative Mariner and
Tokyo Maru geometry remain the next ship-data milestones.

For the canonical ellipsoid with ``(a,b,c)=(3,1,0.8)`` m, the analytical
lower-half sway added mass is ``3340.618`` kg at ``\rho=1000`` kg/m³. The
constant-panel indirect sequence used by `solve_wang_maneuvering` reduces its
error from ``+7.87\%`` at 32 panels to ``+5.84\%`` at 200 panels; the direct
formulation used by `solve_potential_flow_maneuvering` reaches ``-0.39\%`` at
the same 200 panels. An independent 128-panel radiation calculation gives
``3321.805`` kg at ``\omega=0.01`` rad/s, a ``-0.563\%`` analytical error, with
damping at roundoff scale. This completes the analytical and low-frequency
portions of Stage 3. The wall-sided canonical hull, systematic stern-cutoff
sensitivity, and raw captive-test ingestion remain open validation gates.

## Corrections applied after verification

Three defects were identified by the verification work in
`notebooks/01_maneuvering_theory.ipynb`. All three are now addressed. The
previous behaviour remains reachable through keyword arguments so that earlier
results can be reproduced exactly.

### 1. The convective term was an incomplete linearization

Wang's pressure ``p = -\rho(\partial_t - U\partial_x)\Phi`` linearizes the
steady Bernoulli equation about the **undisturbed uniform stream** rather than
about the double-body base flow. Working in the body frame, where the
quasi-steady flow is stationary, the complete linearization for a mode ``j``
with rigid-body velocity field ``\boldsymbol{u}_j`` is

```math
\frac{\partial p}{\partial v_j}
= \rho U\left(\frac{\partial\phi_j}{\partial x}
+\boldsymbol{u}_j\cdot\nabla\phi_1
-\nabla\phi_1\cdot\nabla\phi_j\right),
```

of which only ``\partial\phi_j/\partial x`` was retained. The dropped terms are
first order in ``v_j`` and ``O(U)`` — the same order as the retained one. They
are the ``m``-terms of steady-forward-speed theory. Note that the rotation
field is not a gradient, so the yaw case cannot be written as
``\nabla\Phi_0\cdot\nabla\delta\Phi``; the body-frame form above covers both
modes.

Consequence, measured against the exact Kirchhoff result
``N_v = (A_{11}-A_{22})U`` on a triaxial ellipsoid whose added masses are
analytic: the uniform-stream form converges to ``0.908`` times the exact Munk
moment and does not improve with refinement, while the complete form reaches
``0.9992`` at 800 panels. On KVLCC2 the two differ by roughly 14%.

**Resolved.** `steady_pressure_derivative` implements the expression above, and
`solve_potential_flow_maneuvering` selects between them with `base_flow`,
defaulting to `:double_body`. The extra cost is one unit-surge solve, which the
boundary layer already needs. Pass `base_flow=:uniform_stream` for Wang's
published form. `wang_hydrodynamic_derivatives` is deliberately left on the
uniform-stream form, because its purpose is to reproduce the paper.

### 2. The acceleration derivatives used the less accurate formulation

The maneuvering solvers hardcoded `direct=false`. On identical meshes with the
identical Green function and boundary condition, the sway added mass of the
canonical ellipsoid converges as

| panels | indirect | direct |
|---:|---:|---:|
| 32 | ``+7.87\%`` | ``-2.30\%`` |
| 128 | ``+6.79\%`` | ``-0.56\%`` |
| 512 | ``+4.02\%`` | ``-0.18\%`` |

The direct formulation reaches the indirect formulation's best accuracy with
one sixteenth of the panels, and its zero-frequency value matches the radiation
solver's ``\omega\to0`` limit to five digits.

**Resolved.** The indirect path is still needed for the velocity derivatives,
because `evaluate_indirect_potential_gradient` recovers
``\partial\phi/\partial x`` from the source strengths, but the acceleration
derivatives need only ``\phi`` and are now taken from a direct solve of the
same right-hand sides. That is `acceleration_formulation=:direct`, the default
for `solve_potential_flow_maneuvering`; pass `:indirect` to reproduce earlier
results.

### 3. Wang's second truncation station was not implemented

The paper places the Schmitz cut at the maximum-area section for ``Y_v'`` and
``N_v'`` but explicitly further aft for ``Y_r'`` and ``N_r'``, and recommends
Clarke's regressions

```math
Y_r' = \left(-\tfrac12 + 1.73\tfrac{B}{L}\right)Y_v',
\qquad
N_r' = 0.08\,Y_v' + \tfrac14\left(1 - 0.938\,C_B\tfrac{B}{L}\right)
\left(-\left(\tfrac{T}{L}\right)^2\right)
```

instead. A single mask applied to all four coefficients is the direct cause of
the incorrect sign of the computed ``Y_r'``.

**Resolved.** `clarke_rotational_derivatives` implements both regressions, and
`rotational_velocity_mask` truncates ``Y_r`` and ``N_r`` at a station
independent of the one used for ``Y_v`` and ``N_v``. For KVLCC2 the Clarke
expressions give ``Y_r'=+0.00382`` and ``N_r'=-0.00255`` against model-test
values of ``+0.005395`` and ``-0.003185`` — correct signs and within 30%, from
the sway derivative alone.

### Note on the reported residuals

`boundary_residual` is ``\|K\sigma - g\|/\|g\|`` with ``\sigma = K\backslash g``,
so it is at roundoff by construction, and contracting the reconstructed gradient
with ``n`` reproduces ``K\sigma`` algebraically. Both confirm that the linear
solve succeeded and neither bounds discretization error. The physical checks are
the analytical ellipsoid, the Munk moment, and the acceleration-reciprocity
residual ``|Y_{\dot r} - N_{\dot v}|``.
