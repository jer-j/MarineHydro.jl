using MarineHydro
using PyCall

cpt = pyimport("capytaine")

# Replace this canonical body with the immersed portion of a ship-hull mesh.
length = 6.0
radius = 1.0
cptmesh = cpt.meshes.predefined.mesh_horizontal_cylinder(
    length=length,
    radius=radius,
    center=(0.0, 0.0, 0.0),
    faces_max_radius=0.35,
).immersed_part()
mesh = Mesh(cptmesh)

forward_speed = 2.0
rho = 1025.0

# The default calculation uses the complete wetted hull for both the
# acceleration-dependent and velocity-dependent derivatives.
whole_hull = solve_wang_maneuvering(mesh, forward_speed; rho)

# Wang's Schmitz correction omits the separated stern contribution from the
# velocity-dependent integrals. Here x points toward the bow, so panels aft of
# x = -1.5 m are excluded. The potential only needs to be solved once.
velocity_mask = wang_stern_mask(mesh, -1.5; bow_direction=:positive_x)
truncated = wang_hydrodynamic_derivatives(
    mesh,
    whole_hull.potential,
    whole_hull.potential_gradient[:, 1],
    forward_speed;
    rho,
    velocity_mask,
)

whole_hull_prime = nondimensionalize_wang_derivatives(
    whole_hull.derivatives,
    length,
    forward_speed;
    rho,
)
truncated_prime = nondimensionalize_wang_derivatives(
    truncated,
    length,
    forward_speed;
    rho,
)

println("Relative boundary residual: ", whole_hull.boundary_residual)
println("Whole-hull dimensional derivatives: ", whole_hull.derivatives)
println("Truncated dimensional derivatives: ", truncated)
println("Whole-hull nondimensional derivatives: ", whole_hull_prime)
println("Truncated nondimensional derivatives: ", truncated_prime)
