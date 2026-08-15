# Wang et al. (2000) validation hulls

Wang, Bakountouzis, and Katory validate their method with the Mariner and
Tokyo Maru models. `hulls.toml` records the dimensions, test Froude numbers,
panel counts, and Schmitz integration limits given in the paper.

The paper contains body-plan figures but not numerical station offsets. No
authoritative public offset table has yet been located for either exact model.
Consequently, this repository does not manufacture a hull from an uncalibrated
figure or label a generic Mariner-like form as the experimental geometry.

The remaining reproduction work is:

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
