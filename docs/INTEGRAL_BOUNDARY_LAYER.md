# Integral boundary-layer correction

## Implemented model

MarineHydro now contains a quasi-3D viscous-inviscid interaction model intended
for unrestricted, double-body maneuvering calculations. It preserves the
longitudinal topology of structured Gothenburg hull patches, marches Head's
two-equation turbulent integral boundary layer over many strips on the full 3D
hull surface, and returns both wall shear and displacement transpiration.

For each strip, the attached turbulent layer satisfies

```math
\frac{d\theta}{ds}+(H+2)\frac{\theta}{U_e}\frac{dU_e}{ds}
=\frac{C_f}{2},
\qquad
\frac{d}{ds}(U_e\theta H_1)=U_eE.
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
v_n^{BL}=\frac{1}{w}\frac{d}{ds}(wU_e\delta^*).
```

Published formulations may show the opposite sign when either the body-normal
or streamline-coordinate convention is reversed.

This is deliberately a correction to the whole-hull inviscid result. It should
not be added to a Schmitz-truncated result because the truncation is already a
semi-empirical viscous and vortex correction.

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

For immediate integrated-load validation, the script
`validation/gothenburg2010/kvlcc2_boundary_layer_validation.jl` compares the
predicted linear derivatives with the model-test-derived MMG values reported
for KVLCC2. It reports whole-hull potential flow, Schmitz potential flow, direct
Head shear, Head displacement pressure, and the combined whole-hull result as
separate rows. No experimental derivative is used to tune the closure.

Run:

```sh
julia --project=. \
  validation/gothenburg2010/kvlcc2_boundary_layer_validation.jl
```

Use `--fine` to add a fourth surface resolution.

For a coupled diagnostic at one resolution:

```sh
julia --project=. \
  validation/gothenburg2010/kvlcc2_boundary_layer_validation.jl \
  --shape=8x5 --coupling-iterations=8 \
  --output=validation/gothenburg2010/results/kvlcc2_boundary_layer_coupled.csv
```

### Current result

The current comparison is a validation failure in the scientific sense, not a
calibration. At the finest one-way resolution, the predicted and reference
Wang-normalized derivatives are:

| Coefficient | Head-corrected, 480 panels | KVLCC2 MMG reference |
|---|---:|---:|
| $Y_v'$ | $-0.004461$ | $-0.020475$ |
| $Y_r'$ | $-0.004760$ | $+0.005395$ |
| $N_v'$ | $-0.015075$ | $-0.008905$ |
| $N_r'$ | $-0.001367$ | $-0.003185$ |

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
- Altair, *FlightStream Theory Manual*,
  <https://flightstream-theory.altair.com/>.
