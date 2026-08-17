using BenchmarkTools
using LinearAlgebra
using MarineHydro
using StaticArrays

# Hot paths of the surface boundary-layer solve, on a real hull.
#
# The quantities to watch are allocations rather than time: `cell_residual` is
# the innermost kernel and must stay at zero, and `build_surface_cache` is
# rebuilt on every evaluation of a differentiated residual, so an allocation
# there is paid thousands of times per solve. A regression in either shows up
# here long before it shows up in a wall-clock figure.

const DATA = joinpath(@__DIR__, "..", "validation", "gothenburg2010", "data",
    "KVLCC2")
const PATHS = [joinpath(DATA, "kvlcc_bow1.dat"), joinpath(DATA, "kvlcc2_stn1.dat")]

suite = BenchmarkGroup()

if all(isfile, PATHS)
    forward_speed = 0.142 * sqrt(SETTINGS.g)
    viscosity = forward_speed / 4.6e6
    grid = read_gothenburg2010_panel_grid(PATHS; target_shape = (16, 9))
    mesh = grid.mesh
    surge = solve_rigid_body_potential(mesh, :surge)
    edge_velocity = body_relative_edge_velocity(mesh, forward_speed, 0.0, 0.0;
        surge_gradient = surge.potential_gradient)
    topology = build_surface_topology(mesh)
    metrics = build_surface_metrics(mesh, topology)
    cache = build_surface_cache(mesh, edge_velocity, viscosity; topology, metrics)
    order = flow_ordering(cache)
    inflow = inflow_states(mesh, cache)
    states = initial_states(mesh, cache, order, inflow)
    shear = MarineHydro._lagged_shear(states, mesh, cache, order)
    donors = MarineHydro._donor_states(states, 1, cache, inflow)
    own = SVector{3}(@view states[:, 1])

    suite["geometry"] = BenchmarkGroup()
    suite["geometry"]["build_surface_topology"] = @benchmarkable build_surface_topology($mesh)
    suite["geometry"]["build_surface_metrics"] = @benchmarkable build_surface_metrics($mesh,
        $topology)

    suite["assembly"] = BenchmarkGroup()
    suite["assembly"]["build_surface_cache"] = @benchmarkable build_surface_cache($mesh,
        $edge_velocity, $viscosity; topology = $topology, metrics = $metrics)
    suite["assembly"]["global_residual"] = @benchmarkable global_residual($states,
        $mesh, $cache, $inflow, $shear)
    suite["assembly"]["assemble_jacobian"] = @benchmarkable assemble_jacobian($states,
        $mesh, $cache, $inflow, $shear)

    suite["kernel"] = BenchmarkGroup()
    suite["kernel"]["cell_residual"] = @benchmarkable cell_residual($own, $donors, 1,
        $mesh, $cache, $shear)
    suite["kernel"]["flow_ordering"] = @benchmarkable flow_ordering($cache)

    suite["solve"] = BenchmarkGroup()
    suite["solve"]["surface_boundary_layer"] = @benchmarkable solve_surface_boundary_layer($mesh,
        $edge_velocity, $viscosity; rho = 1025.0, topology = $topology,
        metrics = $metrics)
else
    @info "skipping surface boundary-layer benchmarks; run fetch_geometry.sh first"
end

suite
