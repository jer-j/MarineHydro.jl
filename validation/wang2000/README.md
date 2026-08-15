# Wang et al. (2000) validation hulls

Wang, Bakountouzis, and Katory validate their method with the Mariner and
Tokyo Maru models. `hulls.toml` records the dimensions, test Froude numbers,
panel counts, and Schmitz integration limits given in the paper.

## Canonical analytical validation

`ellipsoid_convergence.jl` supplies the Stage 3 benchmark that does not depend
on unavailable ship offsets. For a triaxial ellipsoid, it evaluates the
classical potential coefficients

```math
\alpha_i=abc\int_0^\infty
\frac{\mathrm{d}\lambda}
{(a_i^2+\lambda)
\sqrt{(a^2+\lambda)(b^2+\lambda)(c^2+\lambda)}}
```

and added masses

```math
m_i=\rho\frac{4\pi abc}{3}\frac{\alpha_i}{2-\alpha_i}.
```

The lower half ellipsoid at the rigid free surface has sway added mass
$m_y/2$, so Wang's convention requires $Y_{\dot v}=-m_y/2$. Run the default
refinement and low-frequency studies with:

```sh
julia --project=. validation/wang2000/ellipsoid_convergence.jl
```

Add `--fine` for the 288- and 512-panel meshes. The driver writes
`ellipsoid_mesh_convergence.csv` and `ellipsoid_low_frequency.csv` under the
ignored `results/` directory. For $(a,b,c)=(3,1,0.8)\ \mathrm{m}$ and
$\rho=1000\ \mathrm{kg\,m^{-3}}$, the analytical lower-half sway added mass is
$3340.618\ \mathrm{kg}$. The default Wang meshes reduce their relative error
from $+7.87\%$ at 32 panels to $+5.84\%$ at 200 panels. On the independent
128-panel radiation calculation, the $\omega=0.01\ \mathrm{rad\,s^{-1}}$
result is $3321.805\ \mathrm{kg}$, or $-0.563\%$ relative to the analytical
value, while the radiation damping tends to zero.

The corresponding plots and executable calculations are documented in
[`notebooks/wang_ellipsoid_validation.ipynb`](../../notebooks/wang_ellipsoid_validation.ipynb).

## Wang ship-model reproduction

The paper contains body-plan figures but not numerical station offsets. No
authoritative public offset table has yet been located for either exact model.
Consequently, this repository does not manufacture a hull from an uncalibrated
figure or label a generic Mariner-like form as the experimental geometry.

The remaining paper-specific reproduction work is:

1. Obtain archival offsets from the original experimental sources, or digitize
   the published body plans with an explicit image-calibration uncertainty.
2. Generate the reported $2\times147$ and $2\times170$ panelizations.
3. Apply the paper's velocity-integral limits: forward perpendicular to
   station 4 for $Y_v'$ and $N_v'$, and to station 2 for $Y_r'$ and $N_r'$.
4. Compare whole-hull, truncated, and independent-yaw results against the
   plotted experimental and theoretical curves.

The primary Mariner experiment is Chislett and Strøm-Tejsen, “Planar Motion
Mechanism Tests and Full-Scale Steering and Manoeuvring Predictions for a
Mariner Class Vessel,” *International Shipbuilding Progress* 12(129), 1965,
[DOI 10.3233/ISP-1965-1212902](https://doi.org/10.3233/ISP-1965-1212902).
