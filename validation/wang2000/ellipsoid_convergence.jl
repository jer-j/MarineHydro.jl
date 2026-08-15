using MarineHydro

const DEFAULT_SEMI_AXES = (3.0, 1.0, 0.8)
const DEFAULT_DENSITY = 1000.0
const DEFAULT_FORWARD_SPEED = 1.0
const DEFAULT_FREQUENCIES = (0.2, 0.1, 0.05, 0.02, 0.01)

"""
    has_flag(flag)

Return `true` when `flag` appears in the command-line arguments.
"""
has_flag(flag) = flag in ARGS

"""
    argument_value(prefix, default)

Return the value following `prefix` in the first matching command-line
argument.
"""
function argument_value(prefix, default)
    argument_index = findfirst(argument -> startswith(argument, prefix), ARGS)
    return isnothing(argument_index) ? default : split(
        ARGS[argument_index],
        "=";
        limit=2,
    )[2]
end

function parse_semi_axes(value)
    values = parse.(Float64, split(value, ','))
    length(values) == 3 || error("--semi-axes must have the form a,b,c")
    return Tuple(values)
end

shape_label(shape) = "$(shape[1])x$(shape[2])"

function write_mesh_convergence(path, semi_axes, density, shapes)
    analytical_added_mass = ellipsoid_added_mass(
        semi_axes;
        density,
    ).y / 2
    analytical_volume = 2pi * prod(semi_axes) / 3
    mkpath(dirname(path))
    open(path, "w") do io
        println(
            io,
            "shape,faces,volume,analytical_volume,volume_relative_error," *
            "wang_added_mass,analytical_added_mass,added_mass_relative_error," *
            "boundary_residual,Y_rdot_over_Y_vdot,reciprocity_error",
        )
        for shape in shapes
            mesh = surface_piercing_ellipsoid_mesh(
                semi_axes;
                longitudinal_panels=shape[1],
                girth_panels=shape[2],
            )
            result = solve_wang_maneuvering(
                mesh,
                DEFAULT_FORWARD_SPEED;
                rho=density,
            )
            volume = mesh_signed_volume(mesh)
            wang_added_mass = -result.derivatives.Y_vdot
            row = (
                shape_label(shape),
                mesh.nfaces,
                volume,
                analytical_volume,
                volume / analytical_volume - 1,
                wang_added_mass,
                analytical_added_mass,
                wang_added_mass / analytical_added_mass - 1,
                result.boundary_residual,
                result.derivatives.Y_rdot / result.derivatives.Y_vdot,
                result.derivatives.N_vdot - result.derivatives.Y_rdot,
            )
            println(io, join(row, ','))
            println(row)
        end
    end
    return nothing
end

function write_low_frequency_limit(
    path,
    semi_axes,
    density,
    shape,
    frequencies,
)
    mesh = surface_piercing_ellipsoid_mesh(
        semi_axes;
        longitudinal_panels=shape[1],
        girth_panels=shape[2],
    )
    zero_frequency = -solve_wang_maneuvering(
        mesh,
        DEFAULT_FORWARD_SPEED;
        rho=density,
    ).derivatives.Y_vdot
    analytical = ellipsoid_added_mass(semi_axes; density).y / 2
    previous_density = SETTINGS.rho
    mkpath(dirname(path))
    try
        set_rho!(density)
        open(path, "w") do io
            println(
                io,
                "shape,faces,omega,radiation_added_mass,radiation_damping," *
                "wang_zero_frequency,analytical_added_mass," *
                "radiation_to_wang_relative_difference," *
                "radiation_to_analytical_relative_error",
            )
            for omega in frequencies
                radiation = calculate_radiation_forces(
                    mesh,
                    [0.0, 1.0, 0.0],
                    omega,
                )
                row = (
                    shape_label(shape),
                    mesh.nfaces,
                    omega,
                    radiation[1],
                    radiation[2],
                    zero_frequency,
                    analytical,
                    radiation[1] / zero_frequency - 1,
                    radiation[1] / analytical - 1,
                )
                println(io, join(row, ','))
                println(row)
            end
        end
    finally
        set_rho!(previous_density)
    end
    return nothing
end

function main()
    fine = has_flag("--fine")
    semi_axes = parse_semi_axes(argument_value(
        "--semi-axes=",
        join(DEFAULT_SEMI_AXES, ','),
    ))
    output_directory = argument_value(
        "--output-dir=",
        joinpath(@__DIR__, "results"),
    )
    shapes = fine ?
        ((8, 4), (12, 6), (16, 8), (20, 10), (24, 12), (32, 16)) :
        ((8, 4), (12, 6), (16, 8), (20, 10))
    radiation_shape = fine ? (24, 12) : (16, 8)

    write_mesh_convergence(
        joinpath(output_directory, "ellipsoid_mesh_convergence.csv"),
        semi_axes,
        DEFAULT_DENSITY,
        shapes,
    )
    write_low_frequency_limit(
        joinpath(output_directory, "ellipsoid_low_frequency.csv"),
        semi_axes,
        DEFAULT_DENSITY,
        radiation_shape,
        DEFAULT_FREQUENCIES,
    )
end

main()
