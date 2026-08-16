# Integral boundary-layer correction

## Implemented model

MarineHydro now contains a quasi-3D viscous-inviscid interaction model intended
for unrestricted, double-body maneuvering calculations. It preserves the
longitudinal topology of structured Gothenburg hull patches, marches Head's
two-equation turbulent integral boundary layer over many strips on the full 3D
hull surface, and returns both wall shear and displacement transpiration.

For each strip, the attached turbulent layer satisfies

```math
\f\frac{d\theta}{ds}+(H+2)\f\frac{\theta}{U_e}\f\frac{dU_e}{ds}
=\f\frac{C_f}{2},
\qquad
\f\frac{d}{ds}(U_e\theta H_1)=U_eE.
```

The edge velocity combines the unit surge, sway, and yaw BEM potentials with
the body-relative rigid velocity. Forward-mode automatic differentiation gives
the derivatives of wall shear and transpiration with respect to sway velocity
and yaw rate while the layer is attached and coupling is one-way. A symmetric
linearization is selected automatically when separation or iterative BEM
coupling introduces a nonsmooth or complex-valued solve path. The
differentiated transpiration is solved through the same BEM operator and
contributes an additional pressure load.

Set `coupling_iterations` to a positive integer to feed the tangential velocity
from the transpiration potential back into the integral boundary layer. The
feedback is under-relaxed with `coupling_relaxation`. Zero iterations retains
the reproducible one-way baseline.

With the body normal directed into the fluid and ``s`` increasing from bow to
stern, MarineHydro uses

```math
v_n^{BL}=\f\frac{1}{w}\f\frac{d}{ds}(wU_e\delta^*).
```

Published formulations may show the opposite sign when either the body-normal
or streamline-coordinate convention is reversed.

This is deliberately a correction to the whole-hull inviscid result. It should
not be added to a Schmitz-truncated result because the truncation is already a
semi-empirical viscous and vortex correction.

Following Wang et al., the correction is integrated over the **attached** region
only. Downstream of separation the attached-flow closure is invalid, so any
pressure or shear it predicts there carries no information. The shear
contribution was already attached-only because the closure zeroes the wall
stress past separation; the displacement-pressure integral now is too.

## Truly three-dimensional model: status

`src/three_dimensional_boundary_layer.jl` implements the crossflow-capable
formulation that the quasi-3D model lacks. It solves the full momentum-loss
tensor with two momentum equations and a kinetic energy equation,

```math
\frac{\partial\theta_{xx}}{\partial x}+\frac{\partial\theta_{xy}}{\partial y}+g_x=0,
\quad
\frac{\partial\theta_{yx}}{\partial x}+\frac{\partial\theta_{yy}}{\partial y}+g_y=0,
\quad
\frac{\partial\theta^*_x}{\partial x}+\frac{\partial\theta^*_y}{\partial y}+g_e=0,
```

in a local Cartesian surface basis, with Drela's shear-stress lag advanced
alongside. The tensor structure and the crossflow parameter
``A_c=e^{1-H_k}\beta_w`` are Mughal's, as used by Lokatt and Eller; the scalar
closures for ``H^*``, ``c_f``, ``U_s``, ``C_{\tau,eq}`` and ``C_D`` are Drela's
incompressible turbulent set. Compressibility, transition and the lateral
curvature equation of Drela's four-equation formulation are out of scope.

Two departures from the published sources are deliberate and are recorded here
because they are not transcription:

- The shear-lag source term as printed in Lokatt and Eller's conservation form
  carries a sign that drives the shear stress *away* from equilibrium. The
  physically correct relaxation is used instead, and because it is linear in
  ``\sqrt{C_\tau}`` it is integrated exactly rather than by explicit Euler.
  That matters: on a hull ``\Delta s/\delta\sim40``, and an explicit update
  oscillates violently enough to destroy the solution.
- Two closure relations in that appendix could not be reconciled with the
  quantities they are said to relate, so the corresponding Drela correlations
  are used, which are the same lineage.

### What is validated

On a structured flat plate with a uniform edge velocity the march reproduces
the turbulent flat-plate laws: momentum thickness within 2% and skin friction
within 0.1% of ``\theta=0.036\,x\,Re_x^{-1/5}`` and
``c_f=0.0592\,Re_x^{-1/5}`` at the downstream end, with the crossflow angle and
every crossflow thickness identically zero. That is the essential check that the
tensor closure and the three equations collapse correctly onto the
two-dimensional limit. The first three or four stations carry a startup
transient from the flat-plate initial estimate.

### What does not yet work, and why

The solver **fails on full hull geometry**, and the reason is structural rather
than a matter of tuning. It is a space-marching scheme along the structured
girth-line strips, and marching is well posed only while the external
streamlines stay within the characteristic cone of the marching direction. On
the KVLCC2 forebody the streamlines wrap around the stem at up to about fifty
degrees to the girth lines, so information genuinely propagates from strip to
strip faster than it propagates downstream, and the march breaks down within a
few stations of the bow. Upwinding the girthwise fluxes on the sign of the
crossflow velocity helps marginally but does not fix the underlying
ill-posedness.

The fix is the one both reference implementations use: a **globally coupled
upwind surface solve** rather than a march. Lokatt and Eller discretise the
whole surface with an upwind-biased finite-volume scheme and solve the sparse
nonlinear system by Newton; Zhang uses a finite-element discretisation with a
fully-simultaneous viscous-inviscid Newton coupling. Either replaces the
station-by-station march with a single global system in the four states over
all panels. That is the next piece of work, and it is a substantial one.

Until then the quasi-3D model below remains the one wired into
`viscous_maneuvering_correction`.

## What "quasi-3D" means

The geometry, external velocity, shear direction, panel area, and moment arm
are three-dimensional. Boundary-layer development on each strip is scalar and
uses ``|\boldsymbol{u}_e|``. Consequently, this version does not yet solve a
crossflow momentum-integral equation. A full surface method requires at least
streamwise and crossflow momentum thicknesses, an entrainment equation,
surface metrics and divergence, attachment-line initialization, and wake
continuation.

The classical Cumpsty-Head-Mager family is the intended next extension. An
important precedent for ship hulls is the viscous-inviscid method of Ikehata,
Nagase, and Maruo, which coupled a Hess-Smith potential calculation to an
integral boundary layer and a stern wake correction. These methods are useful
over attached portions of a hull, but published reviews caution that a pure
integral method does not represent a full-form stern wake containing strong
streamwise vortices.

## FlightStream comparison

The FlightStream theory manual provides a useful modern reference architecture.
Its viscous method is not a fully crossflow-coupled 3D integral boundary layer.
It uses two-dimensional integral equations along on-body streamlines on a 3D
surface, explicitly assuming that crossflow terms are not dominant. That is the
same fidelity class as MarineHydro's present quasi-3D model.

Three FlightStream design choices are directly useful here:

1. Trace the inviscid surface vector field to stagnation and attachment lines,
   rather than assuming fixed longitudinal strips.
2. Iterate the displacement-thickness transpiration with the potential solver.
   MarineHydro now supports this through `coupling_iterations`.
3. Treat separation topology and the wake as separate models. FlightStream uses
   phase-plane feature extraction, axial and cylindrical separation criteria,
   and relaxed vortex-wake strands.

The last point cannot be copied literally for a full-form ship. KVLCC2's stern
wake contains a momentum deficit, a strong bilge vortex, and turbulence-driven
streamwise vorticity. A ship-specific wake continuation is needed in addition
to any vortex filament model. The fixed-strip implementation is therefore a
controlled first stage, not the final stern model.

## KVLCC2 validation split

The Gothenburg 2010 Case 1.1a measurements are used to define the physical
condition and future local-flow targets:

- original KVLCC2, double-body, bare hull;
- ``Re=4.6\times10^6`` and ``Fn=0.142``;
- velocity and crossflow maps at ``x/L_{pp}=0.85,0.9825,1.1``;
- total wake fraction at ``r/R=0.4,0.6,0.8,1.0``;
- velocity profiles and Reynolds stresses at the propeller plane.

The measurements show a strong bilge vortex, a secondary counter-rotating
vortex, hook-shaped iso-wake contours, and strong turbulence anisotropy. Those
are validation targets for the future crossflow and wake model, not phenomena
that the present attached strip march can claim to reproduce.

### Closure and integrated-drag verification

Two checks establish that the layer itself is right before any maneuvering
derivative is asked of it. Both are in
[`notebooks/02_kvlcc2_viscous_correction.ipynb`](../notebooks/02_kvlcc2_viscous_correction.ipynb).

1. **Flat plate.** With constant edge speed, the march must recover
   ``\theta = 0.036\,x\,Re_x^{-1/5}`` and ``C_f = 0.0592\,Re_x^{-1/5}``. It sits
   a few per cent below both and closes on them as ``Re_x`` grows — ``\theta``
   from ``7\%`` to ``3\%`` and ``C_f`` from ``9\%`` to ``3\%`` over
   ``7\times10^5 \le Re_x \le 5\times10^6`` — with ``H`` settling at
   ``1.36``–``1.43``.
2. **KVLCC2 friction drag.** The strip march predicts
   ``C_F = 0.947``–``0.950`` times the ITTC-1957 correlation line at
   ``Re = 4.6\times10^6``, essentially unchanged across a sevenfold change in
   panel count (112 to 836). ITTC-57 is deliberately above a true flat-plate
   line, so a bare attached calculation with no form factor landing slightly
   under it is the expected result.

The layer is therefore validated as an attached-flow integral method. What
follows is a failure of the *maneuvering* extraction, not of the closure.

### Derivative comparison

The same notebook compares the predicted linear derivatives with the
model-test-derived MMG values for KVLCC2, reporting whole-hull potential flow,
Schmitz potential flow, direct Head shear, Head displacement pressure, and the
combined whole-hull result as separate rows. No experimental derivative is used
to tune the closure. A coupled diagnostic is available by passing
`coupling_iterations=8` to `viscous_maneuvering_correction`.

### Current result

The current comparison is a validation failure in the scientific sense, not a
calibration. At the finest one-way resolution, the predicted and reference
Wang-normalized derivatives are:

| Coefficient | Head-corrected, 480 panels | KVLCC2 MMG reference |
|---|---:|---:|
| $Y_v'$ | $-0.004037$ | $-0.020475$ |
| $Y_r'$ | $-0.002999$ | $+0.005395$ |
| $N_v'$ | $-0.012663$ | $-0.008905$ |
| $N_r'$ | $-0.001311$ | $-0.003185$ |

These use the current `base_flow=:double_body` and
`acceleration_formulation=:direct` defaults. Wang's original uniform-stream
linearization gives $N_v'=-0.017252$ on the same mesh against $-0.014839$ for
the double-body form, so roughly $14\,\%$ of the inviscid yaw moment was a
modeling artefact rather than physics.

The 112, 180, 264, and 480-panel sequence is not converged, and the two finest
meshes flag 6 and 10 aft panels as separated. Eight under-relaxed coupling
passes on the 112-panel grid reduce the feedback residual to
$2.39\times10^{-3}$, but move the integrated derivatives farther from the
reference. This isolates the next missing physics: surface-streamline
topology, separation, and momentum-defect and bilge-vortex wake continuation.
It does not support tuning Head's attached-flow closure to the four integrated
coefficients.

## References

- Larsson, Stern, and Visonneau, *Numerical Ship Hydrodynamics*, Gothenburg
  2010 proceedings, <https://doi.org/10.1007/978-94-007-7189-5>.
- Dai and Li, "Manoeuvring Prediction of KVLCC2 with Hydrodynamic Derivatives
  Generated by a Virtual Captive Model Test,"
  <https://doi.org/10.2478/pomr-2019-0062>.
- Ikehata, Nagase, and Maruo, "An Improved Method of Turbulent Boundary Layer
  Theory to Solve Viscous Flow around Ship Stern,"
  <https://doi.org/10.2534/jjasnaoe1968.1982.152_44>.
- Altair, *FlightStream Theory Manual*,
  <https://flightstream-theory.altair.com/>.
- Drela, "XFOIL: An Analysis and Design System,"
  <https://doi.org/10.1007/978-3-642-84010-4_1>.
