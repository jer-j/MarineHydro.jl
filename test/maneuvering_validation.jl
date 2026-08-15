using MarineHydro
using PyCall
using Test

@testset "Wang hemisphere analytical validation" begin
    cpt = pyimport("capytaine")
    radius = 1.0
    rho = 1000.0
    cptmesh = cpt.mesh_sphere(
        radius=radius,
        center=(0, 0, 0),
        resolution=(14, 14),
        name="hemisphere",
    ).immersed_part()
    mesh = Mesh(cptmesh)

    result = solve_wang_maneuvering(mesh, 1.0; rho)

    # Reflection across the rigid free surface completes the hemisphere into a
    # sphere. The sway added mass of a sphere is one half of the displaced-fluid
    # mass, and the force integral on one hemisphere is one half of that value.
    expected_Y_vdot = -rho * pi * radius^3 / 3
    @test result.derivatives.Y_vdot ≈ expected_Y_vdot rtol=0.1
    @test abs(result.derivatives.Y_rdot) < 1e-10 * abs(expected_Y_vdot)
    @test result.derivatives.N_vdot == result.derivatives.Y_rdot
    @test result.boundary_residual < 1e-10
end

@testset "Wang triaxial ellipsoid analytical validation" begin
    semi_axes = (3.0, 1.0, 0.8)
    rho = 1000.0
    mesh = surface_piercing_ellipsoid_mesh(
        semi_axes;
        longitudinal_panels=12,
        girth_panels=6,
    )
    result = solve_wang_maneuvering(mesh, 1.0; rho)

    # Same-sign reflection in the rigid plane completes the lower-half mesh
    # into a full ellipsoid. The wetted-half force integral is therefore half
    # the classical unbounded-fluid sway added mass.
    expected_Y_vdot = -ellipsoid_added_mass(semi_axes; density=rho).y / 2
    @test result.derivatives.Y_vdot ≈ expected_Y_vdot rtol=0.1
    @test abs(result.derivatives.Y_rdot) < 1e-12 * abs(expected_Y_vdot)
    @test result.derivatives.N_vdot == result.derivatives.Y_rdot
    @test result.boundary_residual < 1e-10
end
