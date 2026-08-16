# MarineHydro maneuvering notebooks

These three notebooks are the executable record of the maneuvering-derivative
capability. They replace the earlier standalone drivers under `validation/`:
each notebook now carries the theory, the calculation, the figures and the
interpretation for one stage of the method, and there is no separate script to
keep in sync.

Read them in order.

| Notebook | Contents |
|---|---|
| [`01_maneuvering_theory.ipynb`](01_maneuvering_theory.ipynb) | The zero-frequency double-body boundary-value problem, the indirect boundary-element discretization, Wang's derivative integrals and normalization, the Schmitz truncation, and verification against the analytic triaxial ellipsoid — added mass, the radiation solver's low-frequency limit, the direct-versus-indirect formulation comparison, and an exact ideal-flow check of the velocity derivatives via the Munk moment |
| [`02_kvlcc2_viscous_correction.ipynb`](02_kvlcc2_viscous_correction.ipynb) | Gothenburg 2010 KVLCC2 geometry import and convergence, double-body surface flow with stagnation and streamlines, Head's turbulent integral boundary layer with flat-plate and ITTC-1957 friction-drag verification, and the linearized viscous correction to the velocity derivatives |
| [`03_public_hull_derivatives.ipynb`](03_public_hull_derivatives.ipynb) | KCS, KVLCC2, DTMB 5415, DTC and Wigley under one model: four geometry importers, added masses, static-drift derivatives, and the strict-Wang yaw error as a function of slenderness |
| [`04_three_dimensional_boundary_layer.ipynb`](04_three_dimensional_boundary_layer.ipynb) | The crossflow-capable integral boundary layer: momentum-loss tensor and closure, verification in the two-dimensional limit, local flow field and 3D streamlines on KVLCC2, why a strip march fails on hull geometry, and the Schmitz cut compared with a separation-derived cut |

## Setup

Install the notebook tools once in Julia's shared environment:

```julia
import Pkg
Pkg.activate()
Pkg.add(["IJulia", "Revise"])
```

Then launch Jupyter from the repository root:

```julia
using IJulia
jupyterlab(dir=pwd())
```

The notebooks activate the repository project themselves and load `Revise.jl`
before `MarineHydro`, so source edits are picked up without restarting the
kernel.

## Geometry

Notebooks 2 and 3 need public hull geometry, which is downloaded rather than
committed:

```sh
validation/gothenburg2010/fetch_geometry.sh    # KVLCC2 and KCS (notebooks 2, 3)
validation/public_hulls/fetch_geometry.sh      # adds DTMB 5415 and DTC (notebook 3)
```

Notebook 1 needs no downloads; its geometry is generated analytically.
Notebook 4 needs the Gothenburg geometry.

## Conventions

**Notebooks are committed without outputs.** Every figure is written to the
untracked `results/` directory beside the relevant `validation/` case, and every
downloaded mesh lives in an untracked `data/` directory. Clear outputs before
committing:

```sh
jupyter nbconvert --clear-output --inplace notebooks/*.ipynb
```

Numerical pass/fail criteria live in `test/`. These notebooks are scientific
records — they report what the method does, including where it disagrees with
reference data, and they are not test fixtures.
