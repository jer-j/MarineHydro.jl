import FiniteDifferences
import ForwardDiff
import Zygote

using LinearAlgebra
using MarineHydro
using Test

# Reverse-mode differentiability of the shape-sensitivity chain.
#
# The published solver computes its sensitivities with reverse-mode AD, and that
# is the mode that matters for hull-form work: the design variables are many
# (every vertex, or every shape parameter) and the objective is one number, so a
# single adjoint pass replaces one forward solve per variable. A function that
# only differentiates in forward mode is therefore not finished for this
# purpose, and the usual reason it fails is that it writes into preallocated
# arrays -- Zygote cannot see through `setindex!`.
#
# Each test compares all three of reverse mode, forward mode, and a central
# difference, so a wrong adjoint cannot pass by agreeing with a wrong primal.
@testset "Reverse-mode shape sensitivity" begin
    grid = wigley_hull_panel_grid(3.0, 0.3, 0.1875; longitudinal_panels = 10,
        vertical_panels = 3)
    mesh = grid.mesh
    basis = beam_deformation_basis(mesh)
    central = FiniteDifferences.central_fdm(5, 1)

    function agrees(objective; rtol = 1e-6)
        reverse, = Zygote.gradient(objective, 0.0)
        forward = ForwardDiff.derivative(objective, 0.0)
        difference = central(objective, 0.0)
        @test abs(difference) > 0
        @test reverse≈forward rtol=rtol
        @test reverse≈difference rtol=1e-5
        return reverse
    end

    @testset "deform_mesh" begin
        agrees(amplitude -> mesh_signed_volume(deform_mesh(mesh, basis, amplitude)))
        agrees(amplitude -> mesh_surface_area(deform_mesh(mesh, basis, amplitude)))
    end

    @testset "Sectional geometry and added mass" begin
        agrees() do amplitude
            deformed = StructuredPanelGrid(deform_mesh(mesh, basis, amplitude),
                grid.strips)
            sections = sectional_crossflow_geometry(deformed, falses(mesh.nfaces))
            sum(sectional_added_mass(sections; rho = 1000.0) .* sections.spacing)
        end
    end

    @testset "Shed-vorticity derivative through a deforming hull" begin
        # The quantity the whole exercise exists to differentiate.
        sections0 = sectional_crossflow_geometry(grid, falses(mesh.nfaces))
        station = minimum(sections0.station_x) +
                  0.15 * (maximum(sections0.station_x) - minimum(sections0.station_x))
        gradient = agrees() do amplitude
            deformed = StructuredPanelGrid(deform_mesh(mesh, basis, amplitude),
                grid.strips)
            sections = sectional_crossflow_geometry(deformed, falses(mesh.nfaces))
            shed_vorticity_derivatives(sections, station, 2.0; rho = 1000.0).Y_v
        end
        # Widening the hull does not change the draft, so m22 and hence Y_v are
        # nearly unmoved by a beam deformation; the sensitivity is small but not
        # zero, because the station spacing shifts slightly.
        @test isfinite(gradient)
    end

    @testset "Cross-flow drag through a deforming hull" begin
        agrees() do amplitude
            deformed = StructuredPanelGrid(deform_mesh(mesh, basis, amplitude),
                grid.strips)
            crossflow_drag_loads(deformed, trues(mesh.nfaces), 0.2, 0.0;
                rho = 1000.0).force
        end
    end

    @testset "A draft deformation moves Y_v as the theory says" begin
        # m22 goes as the square of the draft, so scaling the draft by (1 + a)
        # should scale Y_v by (1 + a)^2 and the derivative at zero should be
        # twice the value. This checks the gradient means what it should, not
        # merely that three methods agree on a number.
        draft_basis = zero(mesh.vertices)
        draft_basis[:, 3] .= @view mesh.vertices[:, 3]
        sections0 = sectional_crossflow_geometry(grid, falses(mesh.nfaces))
        station = minimum(sections0.station_x) +
                  0.15 * (maximum(sections0.station_x) - minimum(sections0.station_x))
        objective = function (amplitude)
            deformed = StructuredPanelGrid(deform_mesh(mesh, draft_basis, amplitude),
                grid.strips)
            sections = sectional_crossflow_geometry(deformed, falses(mesh.nfaces))
            shed_vorticity_derivatives(sections, station, 2.0; rho = 1000.0).Y_v
        end
        base = objective(0.0)
        reverse, = Zygote.gradient(objective, 0.0)
        @test reverse≈2 * base rtol=0.05
    end
end
