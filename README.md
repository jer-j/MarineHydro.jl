![Test Status](https://github.com/symbiotic-engineering/MarineHydro.jl/actions/workflows/run_tests.yml/badge.svg)
![GitHub](https://img.shields.io/github/license/symbiotic-engineering/MDOcean)
[![DOI](https://zenodo.org/badge/913583350.svg)](https://doi.org/10.5281/zenodo.19859297)

### Fully Differentiable Boundary Element Solver for Hydrodynamic Sensitivity Analysis of Wave-Structure Interactions

#### Authors

**Kapil Khanal<sup>a,b</sup>, Carlos A. Michelén Ströfer<sup>b</sup>, Matthieu Ancellin<sup>c</sup>, Maha Haji<sup>a</sup>**

#### Affiliations
- **<sup>a</sup>Cornell University**  
  Ithaca, NY 14850, USA
- **<sup>b</sup>Sandia National Laboratories**  
  Albuquerque, NM 87123, USA
- **<sup>c</sup>Eurobios Mews Labs**  
  Paris, France

## 🌟 Research Highlights

- 📚 **Derivation and discussion** of the discrete adjoint method for the boundary integral equations.
- 💻 **Review and implementation** of a differentiable boundary element solver for marine hydrodynamics in Julia.
- 🌊 **Exact gradient calculation** for a pair of floating hemispheres with respect to their dimensions, separation distance, and wave environment.
- ⚡ **Mechanical power optimization** for a pair of wave energy converters using exact gradients.


https://www.sciencedirect.com/science/article/pii/S0141118725002937


Fully-differentiable boundary element solver for marine hydrodynamics. This new solver implements both direct and indirect boundary element formulations and uses two green's function expression, Wu et al, and Delhommeau's varying in their accuracy and speed. 
<img width="632" alt="image" src="https://github.com/user-attachments/assets/16247838-770b-480d-9f2f-d4b0a02054bf" />


> ⚠️ **Note**: This package is **work in progress** 🚧 and a separate public release of the package will be done in the future. This current state of the package contains necessary code to replicate the paper 📄. It will go through a significant change in its API for users in future iterations.
>

### Repository: MarineHydro.jl

#### 📂 Folder Structure

- **📁 .github/workflows**  
  Contains workflow files for automated tasks, such as continuous integration (CI).

- **📊 paper**  
  Includes plots and data generated for the paper.

- **📜 src**  
  Source code files for the `MarineHydro.jl` package, including the main functionality.

- **🧪 test**  
  Contains test files and resources to verify the functionality of the source code.

---

### 🚀 How to Run the Code

1. **Install Julia**  
   Ensure you have Julia installed on your system. You can download it from the [JuliaLang website](https://julialang.org/downloads/).

2. **Clone the Repository**  
   Open a terminal and run:  
   ```bash
   git clone https://github.com/symbiotic-engineering/MarineHydro.jl.git
   cd MarineHydro.jl


3. **Install Dependencies**  
   Start Julia from the terminal in the project directory and run the following:  
   ```julia
   using Pkg
   Pkg.activate(".")
   Pkg.instantiate()
   ```
   You may need to configure PyCall. Note: Running `tests/` automatically install capytaine and uses it. 
4. **Configure PyCall**
   ```
   ENV["PYTHON"] = "/path/to/capytaine_env/bin/python" # run 'which python' in your terminal for this
   ```
  
6. **Using the MarineHydro Module**  
Load the module in your Julia session:  
```julia
using MarineHydro
using PyCall
# import your capytaine mesh
cpt = pyimport("capytaine")
radius = 1.0 #fixed
resolution = (10, 10)
cptmesh = cpt.mesh_sphere(name="sphere", radius=radius, center=(0, 0, 0), resolution=resolution) 
cptmesh.keep_immersed_part(inplace=true)

# declare it Julia mesh
mesh = Mesh(cptmesh)  
ω = 1.03
ζ = [0,0,1] # HEAVE: will be more verbose in future iteration. define it again even if defined in Capytaine.
F = DiffractionForce(mesh,ω,ζ)
A,B = calculate_radiation_forces(mesh,ζ,ω)
```

### Low-frequency maneuvering derivatives

The unrestricted, infinite-depth Wang et al. method is available through
`solve_wang_maneuvering`. It solves a zero-frequency unit-sway problem with the
reflected-Rankine Green function and evaluates the linear sway-yaw acceleration
and velocity derivatives from the surface potential and its longitudinal
gradient.

```julia
result = solve_wang_maneuvering(mesh, forward_speed; rho=1025.0)
derivatives = result.derivatives

velocity_mask = wang_stern_mask(mesh, x_cut)
truncated_derivatives = wang_hydrodynamic_derivatives(
    mesh,
    result.potential,
    result.potential_gradient[:, 1],
    forward_speed;
    rho=1025.0,
    velocity_mask,
)
```

For an independent three-dimensional yaw solution, rather than Wang's
`x * sway_potential` approximation, use:

```julia
result = solve_potential_flow_maneuvering(mesh, forward_speed; rho=1025.0)
sway_potential = result.potentials.sway
yaw_potential = result.potentials.yaw
derivatives = result.derivatives
```

An experimental quasi-3D Head integral boundary-layer correction can be
applied to the whole-hull result. The default is one-way; positive
`coupling_iterations` enable an under-relaxed displacement-transpiration
interaction with the BEM:

```julia
grid = read_gothenburg2010_panel_grid(surface_paths; target_shape=(12, 7))
inviscid = solve_potential_flow_maneuvering(grid.mesh, forward_speed)
surge = solve_rigid_body_potential(grid.mesh, :surge)
viscous = viscous_maneuvering_correction(
    grid,
    surge,
    inviscid,
    forward_speed,
    kinematic_viscosity;
    coupling_iterations=0,
)
corrected = apply_viscous_correction(inviscid.derivatives, viscous)
```

The present boundary-layer state is marched on 3D hull strips but does not yet
solve crossflow momentum or a separated stern wake. See
[`docs/INTEGRAL_BOUNDARY_LAYER.md`](docs/INTEGRAL_BOUNDARY_LAYER.md) for the
equations and validation limits.

Public KVLCC2 and KCS workshop geometry can be fetched and validated with the
scripts in [`validation/gothenburg2010`](validation/gothenburg2010). Historical
Mariner and Tokyo Maru metadata from Wang et al. are recorded in
[`validation/wang2000`](validation/wang2000). The KVLCC2 MMG normalization and
the current automatic Schmitz-cutoff comparison are documented in
[`validation/kvlcc2_maneuvering`](validation/kvlcc2_maneuvering).

See [`examples/wang_open_deep_water.jl`](examples/wang_open_deep_water.jl) for
a complete calculation and
[`docs/WANG_MANEUVERING_PLAN.md`](docs/WANG_MANEUVERING_PLAN.md) for the
governing equations, validation gates, and deferred scope.

6. **Differentiability** :
For differentiability with respect to mesh dimension, use `paper/MeshGradients_singlebody.jl`
Differentiability needs an AD engine: use Zygote
```julia
using Zygote
A_w_grad, = Zygote.gradient(w -> calculate_radiation_forces(mesh,ζ,w)[1],ω)
```
