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
