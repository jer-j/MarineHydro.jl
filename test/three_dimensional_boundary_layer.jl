using LinearAlgebra
using MarineHydro
using StaticArrays
using Test

"""
    flat_plate_grid(; kwargs...)

A structured flat plate. With a uniform edge velocity along the strips the flow
is two-dimensional, so the three-dimensional march must collapse onto the
classical turbulent flat-plate laws.
"""
function flat_plate_grid(; length_x = 5.0, width_y = 1.0, stations = 120, strips = 9)
    x = range(0.05, length_x; length = stations + 1)
    y = range(-width_y / 2, width_y / 2; length = strips + 1)
    vertex_rows = NTuple{3, Float64}[]
    face_rows = NTuple{4, Int}[]
    center_rows = NTuple{3, Float64}[]
    normal_rows = NTuple{3, Float64}[]
    areas = Float64[]
    radii = Float64[]
    strip_lists = [Int[] for _ in 1:strips]
    for j in 1:strips, i in 1:stations
        corners = ((x[i], y[j], 0.0), (x[i + 1], y[j], 0.0),
            (x[i + 1], y[j + 1], 0.0), (x[i], y[j + 1], 0.0))
        base = length(vertex_rows)
        append!(vertex_rows, corners)
        push!(face_rows, (base, base + 1, base + 2, base + 3))
        push!(center_rows, ((x[i] + x[i + 1]) / 2, (y[j] + y[j + 1]) / 2, 0.0))
        push!(normal_rows, (0.0, 0.0, 1.0))
        area = (x[i + 1] - x[i]) * (y[j + 1] - y[j])
        push!(areas, area)
        push!(radii, sqrt(area))
        push!(strip_lists[j], length(face_rows))
    end
    matrix(rows) = permutedims(reduce(hcat, [collect(row) for row in rows]))
    mesh = MarineHydro.Mesh(matrix(vertex_rows), matrix(face_rows),
        matrix(center_rows), matrix(normal_rows), areas, radii,
        length(vertex_rows), length(face_rows))
    return StructuredPanelGrid(mesh, strip_lists)
end

@testset "Three-dimensional integral boundary layer" begin
    @testset "Closure argument validation" begin
        @test_throws ArgumentError ThreeDimensionalClosure(minimum_shape_factor = 0.9)
        @test_throws ArgumentError ThreeDimensionalClosure(initial_shape_factor = 1.0)
        @test_throws ArgumentError ThreeDimensionalClosure(separation_shape_factor = 1.2)
        @test_throws ArgumentError ThreeDimensionalClosure(maximum_crossflow_angle = -1)
        closure = ThreeDimensionalClosure()
        @test closure.initial_shape_factor < closure.separation_shape_factor
    end

    @testset "Scalar closures against known turbulent values" begin
        # A zero-pressure-gradient turbulent layer at Re_theta = 1000 sits near
        # H = 1.4, H* = 1.75, c_f = 0.0042.
        @test turbulent_kinetic_shape_factor(1.4, 1000.0)≈1.75 rtol=0.05
        @test turbulent_skin_friction(1.4, 1000.0)≈0.0042 rtol=0.05
        # Skin friction falls with Reynolds number and with shape factor.
        @test turbulent_skin_friction(1.4, 10_000.0) < turbulent_skin_friction(1.4, 1000.0)
        @test turbulent_skin_friction(2.0, 1000.0) < turbulent_skin_friction(1.4, 1000.0)
        # The kinetic energy shape factor falls as the profile fills out.
        @test turbulent_kinetic_shape_factor(2.0, 1000.0) <
              turbulent_kinetic_shape_factor(1.4, 1000.0)
        slip = turbulent_slip_velocity(1.4, turbulent_kinetic_shape_factor(1.4, 1000.0))
        @test 0 < slip < 1
        @test equilibrium_shear_coefficient(1.4,
            turbulent_kinetic_shape_factor(1.4, 1000.0), slip) > 0
        # The layer thickness is several times the momentum thickness.
        @test layer_thickness(1.0, 1.4) > 5
    end

    @testset "Equilibrium shape factor" begin
        for reynolds in (500.0, 2000.0, 20_000.0)
            shape = equilibrium_shape_factor(reynolds)
            @test 1.2 < shape < 1.8
            kinetic = turbulent_kinetic_shape_factor(shape, reynolds)
            friction = turbulent_skin_friction(shape, reynolds)
            slip = turbulent_slip_velocity(shape, kinetic)
            shear = equilibrium_shear_coefficient(shape, kinetic, slip)
            dissipation = 0.5 * friction * slip + shear * (1 - slip)
            # By construction the energy equation is in balance there.
            @test 2 * dissipation≈kinetic * friction / 2 rtol=1e-6
        end
    end

    @testset "Crossflow closure reduces to two dimensions" begin
        @test crossflow_parameter(1.4, 0.0) == 0
        @test crossflow_parameter(1.4, 0.1) > 0
        @test crossflow_parameter(1.4, -0.1) < 0
        # A larger shape factor weakens the crossflow response.
        @test abs(crossflow_parameter(2.0, 0.1)) < abs(crossflow_parameter(1.4, 0.1))

        thickness, shape, kinetic = 0.002, 1.4, 1.75
        tensors = thickness_tensors(thickness, shape, kinetic, 0.0)
        @test tensors.momentum[1, 1] == thickness
        @test tensors.momentum[1, 2] == 0
        @test tensors.momentum[2, 1] == 0
        @test tensors.momentum[2, 2] == 0
        @test tensors.displacement[1] == shape * thickness
        @test tensors.displacement[2] == 0
        @test tensors.energy[1] == kinetic * thickness
        @test tensors.energy[2] == 0

        # With crossflow every off-diagonal component switches on.
        skewed = thickness_tensors(thickness, shape, kinetic, 0.2)
        @test skewed.momentum[2, 1] != 0
        @test skewed.momentum[1, 2] != 0
        @test skewed.displacement[2] != 0
        @test skewed.energy[2] != 0
        # The streamwise momentum thickness is unaffected by definition.
        @test skewed.momentum[1, 1] == thickness
    end

    @testset "Basis rotation" begin
        tensor = @SMatrix [1.0 2.0; 3.0 4.0]
        @test MarineHydro.rotate_to_local(tensor, 1.0, 0.0) ≈ tensor
        vector = @SVector [1.0, 2.0]
        @test MarineHydro.rotate_vector_to_local(vector, 1.0, 0.0) ≈ vector
        # A quarter turn maps the streamwise axis onto the crossflow axis.
        rotated = MarineHydro.rotate_vector_to_local(@SVector([1.0, 0.0]), 0.0, 1.0)
        @test rotated ≈ [0.0, 1.0]
        # Rotation preserves the trace of a second-rank tensor.
        angle = 0.7
        turned = MarineHydro.rotate_to_local(tensor, cos(angle), sin(angle))
        @test tr(turned)≈tr(tensor) atol=1e-12
    end

    @testset "Shear-stress lag is unconditionally stable" begin
        state = (shape_factor = 1.4, kinetic = 1.75, slip = 0.55,
            momentum_thickness = 1.0e-4)
        equilibrium = equilibrium_shear_coefficient(state.shape_factor, state.kinetic,
            state.slip)
        thickness = layer_thickness(state.momentum_thickness, state.shape_factor)
        # A step far longer than the layer thickness must relax onto
        # equilibrium rather than oscillate, which explicit Euler would not.
        far = advance_shear_coefficient(1e-8, state, 1000 * thickness)
        @test far≈equilibrium rtol=1e-6
        from_above = advance_shear_coefficient(10 * equilibrium, state,
            1000 * thickness)
        @test from_above≈equilibrium rtol=1e-6
        # A short step moves partway and stays bounded by the two states.
        near = advance_shear_coefficient(equilibrium / 4, state, thickness / 10)
        @test equilibrium / 4 < near < equilibrium
    end

    @testset "Marching topology" begin
        grid = flat_plate_grid(stations = 20, strips = 5)
        topology = build_march_topology(grid)
        @test topology.station_count == 20
        @test sort(topology.girth_order) == 1:5
        # The girth ordering must chain neighbours, so consecutive entries
        # differ by one strip index on a plate.
        @test all(abs.(diff(topology.girth_order)) .== 1)
        @test all(topology.streamwise_spacing .> 0)
        @test all(topology.girth_spacing .> 0)
        @test all(topology.arclength[2, :] .> topology.arclength[1, :])
        # The local basis is orthonormal and tangent to the plate.
        for station in 1:20, position in 1:5
            tangent = topology.tangent[station, position, :]
            binormal = topology.binormal[station, position, :]
            @test norm(tangent)≈1 atol=1e-12
            @test norm(binormal)≈1 atol=1e-12
            @test dot(tangent, binormal)≈0 atol=1e-12
        end
    end

    @testset "Flat plate reduces to the two-dimensional laws" begin
        grid = flat_plate_grid()
        mesh = grid.mesh
        viscosity = 1e-6
        edge_speed = 1.0
        edge_velocity = zeros(mesh.nfaces, 3)
        edge_velocity[:, 1] .= edge_speed
        layer = solve_three_dimensional_boundary_layer(grid, edge_velocity, viscosity;
            rho = 1000.0)

        # No crossflow can be generated by a two-dimensional external flow.
        @test maximum(abs, layer.crossflow_angle) == 0
        @test all(iszero, layer.crossflow_displacement_thickness)
        @test !any(layer.separated)
        @test all(layer.attached)

        # Downstream, where the initial transient has washed out, the march
        # must sit on the 1/5-power laws.
        middle_strip = grid.strips[5]
        for station in (80, 120)
            panel = middle_strip[station]
            x = mesh.centers[panel, 1]
            reynolds = edge_speed * x / viscosity
            @test layer.momentum_thickness[panel]≈0.036 * x * reynolds^(-0.2) rtol=0.06
            @test layer.skin_friction_coefficient[panel]≈0.0592 * reynolds^(-0.2) rtol=0.06
            @test 1.25 < layer.shape_factor[panel] < 1.45
        end

        # The layer thickens monotonically everywhere. Skin friction is
        # monotone only after the startup transient: the first station is
        # initialised from a flat-plate estimate that is not an exact solution
        # of the marched equations, and it takes three or four stations for the
        # shape factor to settle onto the equilibrium locus.
        thicknesses = [layer.momentum_thickness[panel] for panel in middle_strip]
        @test all(diff(thicknesses) .> 0)
        frictions = [layer.skin_friction_coefficient[panel] for panel in middle_strip]
        @test all(diff(frictions[5:end]) .< 0)
        shapes = [layer.shape_factor[panel] for panel in middle_strip]
        @test all(1.3 .< shapes[5:end] .< 1.45)

        # The traction is the force the fluid exerts on the wall, so it acts
        # along the imposed edge velocity and carries no lateral component.
        @test layer.force[1] > 0
        @test abs(layer.force[2]) < 1e-12 * abs(layer.force[1])
        # Its magnitude is the integrated local friction.
        expected = sum(0.5 * 1000.0 * edge_speed^2 *
                       layer.skin_friction_coefficient[panel] * mesh.areas[panel]
        for panel in 1:mesh.nfaces)
        @test layer.force[1]≈expected rtol=1e-10
        # Displacement transpiration blows outward from an attached layer.
        @test all(layer.transpiration_velocity .> 0)
    end

    @testset "Argument validation" begin
        grid = flat_plate_grid(stations = 10, strips = 3)
        mesh = grid.mesh
        edge_velocity = zeros(mesh.nfaces, 3)
        edge_velocity[:, 1] .= 1.0
        @test_throws DimensionMismatch solve_three_dimensional_boundary_layer(
            grid, zeros(mesh.nfaces, 2), 1e-6)
        @test_throws ArgumentError solve_three_dimensional_boundary_layer(
            grid, edge_velocity, -1.0)
        @test_throws ArgumentError solve_three_dimensional_boundary_layer(
            grid, edge_velocity, 1e-6; rho = 0.0)
        @test_throws ArgumentError solve_three_dimensional_boundary_layer(
            grid, edge_velocity, 1e-6; station_sweeps = 0)
        @test_throws ArgumentError solve_three_dimensional_boundary_layer(
            grid, edge_velocity, 1e-6; march_substeps = 0)
    end

    @testset "Separation-derived truncation mask" begin
        grid = flat_plate_grid(stations = 10, strips = 3)
        mesh = grid.mesh
        separated = falses(mesh.nfaces)
        # Separate the last three stations of the first strip only.
        for station in 8:10
            separated[grid.strips[1][station]] = true
        end
        per_strip = separation_stern_mask(grid, separated; mode = :per_strip)
        @test count(per_strip) == mesh.nfaces - 3
        @test all(per_strip[grid.strips[2]])
        plane = separation_stern_mask(grid, separated; mode = :first_station)
        # A plane cut removes the same stations from every strip.
        @test count(plane) <= count(per_strip)
        @test separation_stern_mask(grid, separated; mode = :none) == trues(mesh.nfaces)
        @test_throws ArgumentError separation_stern_mask(grid, separated; mode = :bogus)
        @test_throws DimensionMismatch separation_stern_mask(grid, falses(3))
        # A margin moves the cut further forward.
        margined = separation_stern_mask(grid, separated; mode = :per_strip, margin = 2)
        @test count(margined) < count(per_strip)
    end
end
