using LinearAlgebra: norm
using MarineHydro
using Printf

const DEFAULT_SHAPE = (32, 17)

function argument_value(prefix, default)
    argument_index = findfirst(argument -> startswith(argument, prefix), ARGS)
    return isnothing(argument_index) ? default : split(
        ARGS[argument_index],
        "=";
        limit=2,
    )[2]
end

function parse_shape(value)
    dimensions = parse.(Int, split(lowercase(value), 'x'))
    length(dimensions) == 2 || error("--shape must have the form NIxNJ")
    return (dimensions[1], dimensions[2])
end

function workshop_paths(data_directory)
    possible_roots = (joinpath(data_directory, "KVLCC2"), data_directory)
    for root in possible_roots
        bow = joinpath(root, "kvlcc_bow1.dat")
        stern = joinpath(root, "kvlcc2_stn1.dat")
        if isfile(bow) && isfile(stern)
            return [bow, stern]
        end
    end
    error(
        "KVLCC2 surface files were not found. Run validation/gothenburg2010/" *
        "fetch_geometry.sh or pass --data-dir=/path/to/data.",
    )
end

function write_number(io, value)
    isfinite(value) || error("visualization output contains a nonfinite value")
    @printf(io, "%.7g", value)
end

function write_number_vector(io, values)
    print(io, '[')
    for (index, value) in enumerate(values)
        index == 1 || print(io, ',')
        write_number(io, value)
    end
    print(io, ']')
end

function write_integer_vector(io, values)
    print(io, '[')
    for (index, value) in enumerate(values)
        index == 1 || print(io, ',')
        print(io, value)
    end
    print(io, ']')
end

function write_streamline(io, line, is_last)
    point_indices = collect(1:2:size(line.points, 1))
    last(point_indices) == size(line.points, 1) ||
        push!(point_indices, size(line.points, 1))
    print(io, "{\"x\":")
    write_number_vector(io, line.points[point_indices, 1])
    print(io, ",\"y\":")
    write_number_vector(io, line.points[point_indices, 2])
    print(io, ",\"z\":")
    write_number_vector(io, line.points[point_indices, 3])
    print(io, ",\"speed\":")
    write_number_vector(io, line.speed[point_indices])
    print(io, '}')
    is_last || print(io, ',')
end

function write_visualization_data(
    output_path,
    shape,
    grid,
    surge,
    streamlines,
    stagnation,
)
    mesh = grid.mesh
    triangle_i = Int[]
    triangle_j = Int[]
    triangle_k = Int[]
    for panel in 1:mesh.nfaces
        a, b, c, d = mesh.faces[panel, :]
        append!(triangle_i, (a - 1, a - 1))
        append!(triangle_j, (b - 1, c - 1))
        append!(triangle_k, (c - 1, d - 1))
    end
    mkpath(dirname(output_path))
    open(output_path, "w") do io
        print(io, "{\"meta\":{\"shape\":[")
        print(io, shape[1], ',', shape[2])
        print(io, "],\"faces\":", mesh.nfaces)
        print(io, ",\"strips\":", length(grid.strips))
        print(io, ",\"boundary_residual\":")
        write_number(io, surge.boundary_residual)
        print(io, ",\"x_positive\":\"bow\",\"incoming_flow\":\"negative_x\"},\"mesh\":{\"x\":")
        write_number_vector(io, mesh.vertices[:, 1])
        print(io, ",\"y\":")
        write_number_vector(io, mesh.vertices[:, 2])
        print(io, ",\"z\":")
        write_number_vector(io, mesh.vertices[:, 3])
        print(io, ",\"i\":")
        write_integer_vector(io, triangle_i)
        print(io, ",\"j\":")
        write_integer_vector(io, triangle_j)
        print(io, ",\"k\":")
        write_integer_vector(io, triangle_k)
        print(io, "},\"streamlines\":[")
        for (index, line) in enumerate(streamlines)
            write_streamline(io, line, index == length(streamlines))
        end
        print(io, "],\"stagnation\":{\"panel\":", stagnation.panel)
        print(io, ",\"point\":")
        write_number_vector(io, stagnation.point)
        print(io, ",\"velocity\":")
        write_number_vector(io, stagnation.velocity)
        print(io, ",\"speed\":")
        write_number(io, stagnation.speed)
        print(io, "}}")
    end
end

function main()
    data_directory = argument_value(
        "--data-dir=",
        joinpath(@__DIR__, "data"),
    )
    shape = parse_shape(argument_value(
        "--shape=",
        "$(DEFAULT_SHAPE[1])x$(DEFAULT_SHAPE[2])",
    ))
    output_path = argument_value(
        "--output=",
        joinpath(
            @__DIR__,
            "results",
            "kvlcc2_surface_flow_$(shape[1])x$(shape[2]).json",
        ),
    )
    paths = workshop_paths(data_directory)
    grid = read_gothenburg2010_panel_grid(paths; target_shape=shape)
    mesh = grid.mesh
    surge = solve_rigid_body_potential(mesh, :surge)
    edge_velocity = body_relative_edge_velocity(
        mesh,
        1.0,
        0.0,
        0.0;
        surge_gradient=surge.potential_gradient,
    )

    longitudinal_extent = maximum(mesh.centers[:, 1]) -
        minimum(mesh.centers[:, 1])
    bow_limit = maximum(mesh.centers[:, 1]) - 0.12 * longitudinal_extent
    bow_mask = BitVector(mesh.centers[:, 1] .>= bow_limit)
    stagnation = surface_stagnation_panel(
        mesh,
        edge_velocity;
        panel_mask=bow_mask,
    )
    stagnation.point[1] > 0 || error(
        "the minimum-speed bow panel is not on the positive-x forebody",
    )
    forebody_mask = mesh.centers[:, 1] .> 0
    sorted_forebody_x_velocity = sort(edge_velocity[forebody_mask, 1])
    median_forebody_x_velocity = sorted_forebody_x_velocity[
        cld(length(sorted_forebody_x_velocity), 2)
    ]
    median_forebody_x_velocity < 0 || error(
        "the complete body-relative velocity does not point toward negative x",
    )

    seed_panels = [first(strip) for strip in grid.strips]
    step_size = 0.003 * longitudinal_extent
    streamlines = trace_surface_streamlines(
        mesh,
        edge_velocity,
        seed_panels;
        step_size,
        max_steps=650,
        neighbor_count=10,
        maximum_surface_distance=0.06 * longitudinal_extent,
    )
    write_visualization_data(
        output_path,
        shape,
        grid,
        surge,
        streamlines,
        stagnation,
    )
    println((
        shape,
        faces=mesh.nfaces,
        strips=length(grid.strips),
        stagnation_point=stagnation.point,
        stagnation_speed_ratio=stagnation.speed,
        surge.boundary_residual,
        output_path,
    ))
end

main()
