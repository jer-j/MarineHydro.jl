module MarineHydro

using DimensionalData
import Graphs
using ForwardDiff
using LinearAlgebra
using LinearAlgebra: cross, dot, norm
using StaticArrays
using ImplicitAD: implicit_linear
import ImplicitAD

const τ̅ = 2π

include("constants.jl")
export SETTINGS, set_g!, set_rho!

# ---------------------------------------------------------------------------
# Green's functions
# ---------------------------------------------------------------------------

include("green_functions/abstract_greens_function.jl")
export greens, gradient_greens, integral, integral_gradient, with_reduced_coordinates
include("green_functions/rankine.jl")
export Rankine
include("green_functions/rankine_reflected.jl")
export RankineReflected
include("green_functions/rankine_reflected_negative.jl")
export RankineReflectedNegative
include("green_functions/wu.jl")
export GFWu
include("green_functions/exact_Guevel_Delhommeau.jl")
export ExactGuevelDelhommeau

# ---------------------------------------------------------------------------
# Geometry: meshes, the surfaces built from them, and their deformation
# ---------------------------------------------------------------------------

include("geometry/meshes.jl")
export Mesh, element, combine_meshes, +, wavebot_mesh

include("geometry/benchmark_meshes.jl")
export StructuredPanelGrid
export read_tecplot_structured_surfaces, read_gothenburg2010_mesh
export read_plot3d_structured_surfaces, structured_surface_waterline_mesh
export read_ascii_stl_waterline_mesh
export read_gothenburg2010_panel_grid
export structured_sectional_area_curve, gothenburg_maximum_section
export mesh_surface_area, mesh_signed_volume
export ellipsoid_potential_coefficients, ellipsoid_added_mass
export surface_piercing_ellipsoid_mesh, surface_piercing_ellipsoid_panel_grid
export wigley_hull_mesh, wigley_hull_panel_grid

include("geometry/mesh_deformation.jl")
export deform_mesh, beam_deformation_basis

include("geometry/surface_topology.jl")
export SurfaceTopology, SurfaceMetrics
export weld_vertices, build_surface_topology, build_surface_metrics
export boundary_edges, interior_edges, connected_components
export panel_adjacency_graph
export edge_normal_for, edge_area_vector_for

# ---------------------------------------------------------------------------
# Boundary element method
# ---------------------------------------------------------------------------

include("bem/bodies.jl")
export FloatingBody, combine_floatingbodies, +

include("bem/problems_and_results.jl")
export LinearPotentialFlowProblem, DiffractionProblem, RadiationProblem
export LinearPotentialFlowResult, DiffractionResult, RadiationResult
export make_result, problems_from_data, assemble_hydrodynamic_coefficients
export create_DimStack, compute_hydrodynamic_coefficients
export compute_and_label_hydrodynamic_coefficients

include("bem/matrix_assembly.jl")
export assemble_matrices, assemble_matrix_wu, solve
export assemble_matrix_ExactGuevelDelhommeau

include("bem/waves.jl")
export FroudeKrylovForce, AiryBC, airy_waves_pressure, airy_waves_velocity
export airy_waves_potential
export radiation_bc, integrate_pressure, compute_bc, compute_wavenumber
export compute_encountered_values
export calculate_radiation_forces, DiffractionForce, diffraction_force

include("bem/solve.jl")
export solve_problem, solve_all_problems

# ---------------------------------------------------------------------------
# Maneuvering: the potential-flow derivatives the boundary layer corrects
# ---------------------------------------------------------------------------

include("maneuvering/maneuvering.jl")
export WangHydrodynamicDerivatives, WangManeuveringResult
export MMGLinearHullDerivatives, mmg_to_wang_velocity_derivatives
export PotentialFlowManeuveringResult
export RigidBodyPotentialResult, solve_rigid_body_potential
export WangRestrictedWaterResult, wang_restricted_water_elevation
export maneuvering_boundary_condition, wang_sway_boundary_condition
export evaluate_indirect_potential_gradient
export steady_pressure_derivative, clarke_rotational_derivatives
export wang_hydrodynamic_derivatives, nondimensionalize_wang_derivatives
export nondimensionalize_maneuvering_derivatives
export potential_flow_hydrodynamic_derivatives
export solve_wang_maneuvering, solve_potential_flow_maneuvering, wang_stern_mask

# ---------------------------------------------------------------------------
# Boundary layer. The surface solver is included before `boundary_layers.jl`
# because `viscous_maneuvering_correction` dispatches on the closure type, so
# `ThreeDimensionalClosure` must exist when its signature is evaluated.
# ---------------------------------------------------------------------------

include("boundary_layer/three_dimensional_boundary_layer.jl")
export ThreeDimensionalClosure, ThreeDimensionalBoundaryLayerResult
export SurfaceMarchTopology, build_march_topology, girth_ordering
export solve_three_dimensional_boundary_layer
export crossflow_parameter, thickness_tensors
export turbulent_kinetic_shape_factor, turbulent_skin_friction
export turbulent_slip_velocity, equilibrium_shear_coefficient, layer_thickness
export equilibrium_shape_factor, advance_shear_coefficient

include("boundary_layer/surface_boundary_layer.jl")
export SurfaceBoundaryLayerCache, SurfaceBoundaryLayerSolution
export build_surface_cache, solve_surface_boundary_layer
export flow_ordering, upwind_graph, upwind_cycles
export inflow_states, initial_states
export cell_residual, global_residual, assemble_jacobian
export attached_flow_domain, surface_curvature
export rotation_into, edge_partner, edge_flux_contribution

include("boundary_layer/boundary_layers.jl")
export HeadTurbulentClosure, IntegralBoundaryLayerResult
export Quasi3DBoundaryLayerResult, ViscousManeuveringDerivatives
export ViscousManeuveringCorrectionResult
export head_kinetic_shape_factor, head_skin_friction_coefficient
export head_entrainment_coefficient, solve_head_boundary_layer
export body_relative_edge_velocity, solve_quasi3d_boundary_layer
export viscous_maneuvering_correction, apply_viscous_correction
export attachment_weight, separation_stern_mask
export nondimensionalize_viscous_derivatives

include("boundary_layer/simultaneous_coupling.jl")
export SimultaneousCouplingResult, solve_simultaneous_coupling
export indirect_gradient_operator

include("boundary_layer/surface_streamlines.jl")
export SurfaceStreamline, surface_stagnation_panel
export trace_surface_streamline, trace_surface_streamlines

# ---------------------------------------------------------------------------
# Shed vorticity and the empirical methods it is judged against
# ---------------------------------------------------------------------------

include("maneuvering/crossflow_drag.jl")
export SectionalCrossflowDrag, sectional_crossflow_geometry, crossflow_drag_loads

include("maneuvering/shed_vorticity.jl")
export sectional_added_mass, shed_vorticity_derivatives, separation_station

include("maneuvering/empirical_derivatives.jl")
export HullParticulars, empirical_sway_derivatives, empirical_methods
export clarke_full_derivatives, hirano_takashina_derivatives

include("maneuvering/virtual_captive_tests.jl")
export VirtualCaptiveRun, VirtualCaptiveTest
export virtual_captive_test, fit_captive_coefficients

end
