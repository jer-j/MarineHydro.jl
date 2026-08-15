# MarineHydro validation notebooks

These notebooks document visualization and validation workflows that are too
expensive or presentation-oriented for the automated test suite. They use a
native Julia kernel, activate the repository project themselves, and load
`Revise.jl` before `MarineHydro.jl` so edits made during an interactive session
are picked up without restarting the kernel.

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

Open one of the following notebooks with the registered Julia kernel:

- `kvlcc2_surface_flow.ipynb` solves the 1,984-panel KVLCC2 surge problem and
  shows the coordinate orientation, bow stagnation region, bow and stern
  surface streamlines without rendering a free-surface plane or elevation.
- `wang_ellipsoid_validation.ipynb` compares Wang's zero-frequency sway added
  mass with the analytical ellipsoid result and the low-frequency limit of
  MarineHydro's radiation solver.
- `public_hull_derivatives.ipynb` solves KCS, KVLCC2, DTMB 5415, DTC, and
  Wigley, renders every waterline-clipped hull, and compares static-drift
  derivatives and sway-yaw added masses with the available public references.

The KVLCC2 notebook expects the checksum-verified Gothenburg geometry:

```sh
validation/gothenburg2010/fetch_geometry.sh
```

The multi-hull notebook fetches all five geometry sources with:

```sh
validation/public_hulls/fetch_geometry.sh
```

Generated figures are written under the ignored validation `results/`
directories. Numerical pass/fail criteria remain in `test/`; notebooks are
executable scientific records rather than test fixtures.
