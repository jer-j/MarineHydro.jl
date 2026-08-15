# Wang et al. (2000) validation hulls

Wang, Bakountouzis, and Katory validate their method with the Mariner and
Tokyo Maru models. `hulls.toml` records the dimensions, test Froude numbers,
panel counts, and Schmitz integration limits given in the paper.

## Canonical analytical validation

[`notebooks/01_maneuvering_theory.ipynb`](../../notebooks/01_maneuvering_theory.ipynb)
supplies the Stage 3 benchmark that does not depend on unavailable ship
offsets. For a triaxial ellipsoid, it evaluates the classical potential
coefficients

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

The lower half ellipsoid at the rigid free surface is exactly the full
ellipsoid in unbounded fluid, so every surface integral is one half of the
corresponding full-body value and Wang's convention requires
$Y_{\dot v}=-m_y/2$. Figures and CSVs are written under the ignored `results/`
directory.

For $(a,b,c)=(3,1,0.8)\ \mathrm{m}$ and $\rho=1000\ \mathrm{kg\,m^{-3}}$ the
analytical lower-half sway added mass is $3340.618\ \mathrm{kg}$. The notebook
establishes four things about it:

1. The indirect solve used by `solve_wang_maneuvering` reduces its relative
   error from $+7.87\%$ at 32 panels only to $+4.02\%$ at 512 panels.
2. The **direct** formulation, on identical meshes with the identical Green
   function, reaches $-0.18\%$ at 512 panels and $-2.30\%$ at 32. The
   acceleration derivatives therefore carry several percent of avoidable
   formulation error.
3. The independent radiation solver at $\omega=0.01\ \mathrm{rad\,s^{-1}}$ gives
   $-0.563\%$ on the 128-panel mesh with damping at roundoff, matching the
   zero-frequency direct solve to five digits.
4. The **velocity** derivatives can be checked exactly, because ideal flow
   requires $Y_v=0$ and $N_v=(A_{11}-A_{22})U$ (the Munk moment). $Y_v$ is zero
   to roundoff. $N_v$ converges to $0.908$ times the exact value and stays
   there under refinement — a modeling error, not a discretization error,
   traced to Wang's incomplete linearization of the steady Bernoulli equation.
   Restoring the dropped $\phi_{1,y}$ and $\nabla\phi_1\!\cdot\!\nabla\phi_v$
   terms recovers the Munk moment to $0.08\%$ at 800 panels.

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
