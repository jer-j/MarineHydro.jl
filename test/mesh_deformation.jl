import FiniteDifferences
import ForwardDiff

using LinearAlgebra
using MarineHydro
using Test

@testset "Mesh deformation and shape sensitivity" begin
    @testset "A zero displacement reproduces the mesh" begin
        grid = wigley_hull_panel_grid(3.0, 0.3, 0.1875; longitudinal_panels = 10,
            vertical_panels = 3)
        mesh = grid.mesh
        same = deform_mesh(mesh, zero(mesh.vertices))
        @test same.vertices ≈ mesh.vertices
        @test same.centers≈mesh.centers rtol=1e-12
        @test same.areas≈mesh.areas rtol=1e-12
        @test same.radii≈mesh.radii rtol=1e-12
        # Normals are recomputed, so they may differ in sign convention only if
        # the importer chose differently; on this hull they must match outright.
        @test same.normals≈mesh.normals rtol=1e-12
        @test_throws DimensionMismatch deform_mesh(mesh, zeros(3, 3))
    end

    @testset "Deforming the beam matches a hull built at that beam" begin
        # The deformation route and the parametric route must agree, which is
        # what says `deform_mesh` recomputes the geometry rather than merely
        # moving points.
        base = wigley_hull_panel_grid(3.0, 0.3, 0.1875; longitudinal_panels = 12,
            vertical_panels = 4).mesh
        amplitude = 0.25
        widened = deform_mesh(base, beam_deformation_basis(base), amplitude)
        direct = wigley_hull_panel_grid(3.0, 0.3 * (1 + amplitude), 0.1875;
            longitudinal_panels = 12, vertical_panels = 4).mesh
        @test widened.vertices≈direct.vertices rtol=1e-12
        @test widened.areas≈direct.areas rtol=1e-12
        @test mesh_signed_volume(widened)≈mesh_signed_volume(direct) rtol=1e-12
        # Only approximately linear in the beam: each panel is a warped
        # quadrilateral whose volume contribution comes from a two-triangle
        # split, and that split's error is itself a function of the aspect
        # ratio the scaling changes.
        @test mesh_signed_volume(widened)≈(1 + amplitude) * mesh_signed_volume(base) rtol=0.01
    end

    @testset "Shape derivative agrees with finite differences" begin
        base = wigley_hull_panel_grid(3.0, 0.3, 0.1875; longitudinal_panels = 12,
            vertical_panels = 4).mesh
        basis = beam_deformation_basis(base)
        volume(amplitude) = mesh_signed_volume(deform_mesh(base, basis, amplitude))
        area(amplitude) = mesh_surface_area(deform_mesh(base, basis, amplitude))
        for quantity in (volume, area)
            forward = ForwardDiff.derivative(quantity, 0.0)
            central = FiniteDifferences.central_fdm(5, 1)(quantity, 0.0)
            @test forward≈central rtol=1e-6
            @test abs(forward) > 0
        end
        # Volume is very nearly linear in a beam scaling, so the derivative sits
        # close to the volume itself; the gap is the two-triangle split error.
        @test ForwardDiff.derivative(volume, 0.0)≈mesh_signed_volume(base) rtol=0.01
    end

    @testset "An imported hull carries dual numbers" begin
        # This is the point of the exercise: the readers build Float64
        # coordinates, so deforming afterwards is what puts a dual-valued mesh
        # on a hull that came out of a file.
        paths = [joinpath(@__DIR__, "..", "validation", "gothenburg2010", "data",
                     "KVLCC2", "kvlcc_bow1.dat"),
            joinpath(@__DIR__, "..", "validation", "gothenburg2010", "data", "KVLCC2",
                "kvlcc2_stn1.dat")]
        if !all(isfile, paths)
            @info "skipping imported-hull deformation test; run fetch_geometry.sh"
        else
            imported = read_gothenburg2010_panel_grid(paths;
                target_shape = (10, 5)).mesh
            @test eltype(imported.vertices) === Float64
            basis = beam_deformation_basis(imported)
            volume(amplitude) = mesh_signed_volume(deform_mesh(imported, basis,
                amplitude))
            forward = ForwardDiff.derivative(volume, 0.0)
            central = FiniteDifferences.central_fdm(5, 1)(volume, 0.0)
            @test forward≈central rtol=1e-6
            @test eltype(deform_mesh(imported, basis,
                ForwardDiff.Dual(0.0, 1.0)).centers) <: ForwardDiff.Dual
        end
    end
end
