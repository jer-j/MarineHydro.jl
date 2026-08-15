using MarineHydro
using TOML

const KVLCC2_LENGTH = 320.0
const KVLCC2_DRAFT = 20.8
const KVLCC2_DISPLACEMENT = 312622.0
const KVLCC2_WETTED_AREA = 27194.0
const KVLCC2_FROUDE_NUMBER = 0.142
const WATER_DENSITY = 1025.0
const MMG_SWAY_ADDED_MASS_PRIME = 0.223
const MMG_YAW_ADDED_INERTIA_PRIME = 0.011
const KVLCC2_MANEUVER_REFERENCE = joinpath(
    @__DIR__,
    "..",
    "kvlcc2_maneuvering",
    "reference.toml",
)

const DERIVATIVE_NAMES = (
    :Y_vdot,
    :Y_v,
    :Y_rdot,
    :Y_r,
    :N_vdot,
    :N_v,
    :N_rdot,
    :N_r,
)

"""
    argument_value(prefix, default)

Return the value from the first command-line argument beginning with `prefix`.
"""
function argument_value(prefix, default)
    argument_index = findfirst(argument -> startswith(argument, prefix), ARGS)
    return isnothing(argument_index) ? default : split(ARGS[argument_index], "="; limit=2)[2]
end

function workshop_paths(data_directory)
    possible_roots = (
        joinpath(data_directory, "KVLCC2"),
        data_directory,
    )
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

function shape_label(shape)
    return isnothing(shape) ? "full" : "$(shape[1])x$(shape[2])"
end

function write_geometry_convergence(path, source_paths)
    reference_area = KVLCC2_WETTED_AREA / KVLCC2_LENGTH^2
    reference_volume = KVLCC2_DISPLACEMENT / KVLCC2_LENGTH^3
    shapes = ((8, 5), (10, 6), (20, 12), (40, 24), nothing)

    open(path, "w") do io
        println(io, "shape,faces,S_over_L2,S_relative_error,volume_over_L3,volume_relative_error")
        for shape in shapes
            mesh = read_gothenburg2010_mesh(source_paths; target_shape=shape)
            area = mesh_surface_area(mesh)
            volume = mesh_signed_volume(mesh)
            row = (
                shape_label(shape),
                mesh.nfaces,
                area,
                (area - reference_area) / reference_area,
                volume,
                (volume - reference_volume) / reference_volume,
            )
            println(io, join(row, ','))
            println(row)
        end
    end
end

function derivative_columns(prefix, derivatives)
    return ["$(prefix)_$(name)" => getproperty(derivatives, name) for name in DERIVATIVE_NAMES]
end

function mmg_velocity_reference()
    reference = TOML.parsefile(KVLCC2_MANEUVER_REFERENCE)
    coefficients = reference["mmg"]["linear_hull"]
    mmg = MMGLinearHullDerivatives(
        coefficients["Y_v"],
        coefficients["Y_r"],
        coefficients["N_v"],
        coefficients["N_r"],
    )
    return mmg_to_wang_velocity_derivatives(
        mmg,
        KVLCC2_DRAFT,
        KVLCC2_LENGTH,
    )
end

function write_velocity_comparison_rows(
    io,
    shape,
    mesh,
    section,
    velocity_mask,
    methods,
    reference,
)
    for (method_name, derivatives) in pairs(methods)
        for coefficient in (:Y_v, :Y_r, :N_v, :N_r)
            bem_value = getproperty(derivatives, coefficient)
            reference_value = getproperty(reference, coefficient)
            difference = bem_value - reference_value
            row = (
                shape_label(shape),
                mesh.nfaces,
                section.x,
                section.maximum_x,
                section.relative_area_tolerance,
                count(velocity_mask),
                method_name,
                coefficient,
                bem_value,
                reference_value,
                difference,
                abs(difference) / abs(reference_value),
            )
            println(io, join(row, ','))
        end
    end
end

function write_derivative_convergence(
    path,
    comparison_path,
    source_paths;
    fine=false,
)
    shapes = fine ? ((8, 5), (10, 6), (12, 7), (16, 9)) : ((8, 5), (10, 6))
    forward_speed = KVLCC2_FROUDE_NUMBER * sqrt(SETTINGS.g)
    draft_ratio = KVLCC2_DRAFT / KVLCC2_LENGTH
    mmg_Y_vdot_reference = -draft_ratio * MMG_SWAY_ADDED_MASS_PRIME
    mmg_N_rdot_reference = -draft_ratio * MMG_YAW_ADDED_INERTIA_PRIME
    mmg_velocity = mmg_velocity_reference()
    maximum_section = gothenburg_maximum_section(source_paths)
    header = String[
        "shape",
        "faces",
        "strict_boundary_residual",
        "exact_sway_boundary_residual",
        "exact_yaw_boundary_residual",
        "exact_acceleration_reciprocity_error",
        "mmg_estimate_Y_vdot",
        "strict_Y_vdot_relative_difference_from_mmg",
        "exact_Y_vdot_relative_difference_from_mmg",
        "mmg_estimate_N_rdot",
        "strict_N_rdot_relative_difference_from_mmg",
        "exact_N_rdot_relative_difference_from_mmg",
    ]
    append!(header, ["strict_$(name)" for name in DERIVATIVE_NAMES])
    append!(header, ["exact_$(name)" for name in DERIVATIVE_NAMES])

    comparison_header = (
        "shape",
        "faces",
        "schmitz_x_cut_over_L",
        "discrete_maximum_x_over_L",
        "relative_area_tolerance",
        "velocity_faces",
        "method",
        "coefficient",
        "bem_wang_prime",
        "mmg_wang_prime",
        "signed_difference",
        "absolute_relative_difference",
    )

    open(path, "w") do io
        open(comparison_path, "w") do comparison_io
            println(io, join(header, ','))
            println(comparison_io, join(comparison_header, ','))
            for shape in shapes
                mesh = read_gothenburg2010_mesh(source_paths; target_shape=shape)
                strict = solve_wang_maneuvering(
                    mesh,
                    forward_speed;
                    rho=WATER_DENSITY,
                )
                exact = solve_potential_flow_maneuvering(
                    mesh,
                    forward_speed;
                    rho=WATER_DENSITY,
                )
                strict_prime = nondimensionalize_maneuvering_derivatives(
                    strict.derivatives,
                    1.0,
                    forward_speed;
                    rho=WATER_DENSITY,
                )
                exact_prime = nondimensionalize_maneuvering_derivatives(
                    exact.derivatives,
                    1.0,
                    forward_speed;
                    rho=WATER_DENSITY,
                )
                velocity_mask = wang_stern_mask(mesh, maximum_section.x)
                strict_schmitz = wang_hydrodynamic_derivatives(
                    mesh,
                    strict.potential,
                    strict.potential_gradient[:, 1],
                    forward_speed;
                    rho=WATER_DENSITY,
                    velocity_mask,
                )
                exact_schmitz = potential_flow_hydrodynamic_derivatives(
                    mesh,
                    exact.potentials.sway,
                    exact.potentials.yaw,
                    exact.potential_gradients.sway[:, 1],
                    exact.potential_gradients.yaw[:, 1],
                    forward_speed;
                    rho=WATER_DENSITY,
                    velocity_mask,
                )
                strict_schmitz_prime = nondimensionalize_maneuvering_derivatives(
                    strict_schmitz,
                    1.0,
                    forward_speed;
                    rho=WATER_DENSITY,
                )
                exact_schmitz_prime = nondimensionalize_maneuvering_derivatives(
                    exact_schmitz,
                    1.0,
                    forward_speed;
                    rho=WATER_DENSITY,
                )
                reciprocity_scale = max(
                    abs(exact.derivatives.Y_rdot),
                    abs(exact.derivatives.N_vdot),
                    eps(Float64),
                )
                reciprocity_error = abs(
                    exact.derivatives.Y_rdot - exact.derivatives.N_vdot,
                ) / reciprocity_scale
                row = Any[
                    shape_label(shape),
                    mesh.nfaces,
                    strict.boundary_residual,
                    exact.boundary_residuals.sway,
                    exact.boundary_residuals.yaw,
                    reciprocity_error,
                    mmg_Y_vdot_reference,
                    abs(
                        (strict_prime.Y_vdot - mmg_Y_vdot_reference) /
                        mmg_Y_vdot_reference,
                    ),
                    abs(
                        (exact_prime.Y_vdot - mmg_Y_vdot_reference) /
                        mmg_Y_vdot_reference,
                    ),
                    mmg_N_rdot_reference,
                    abs(
                        (strict_prime.N_rdot - mmg_N_rdot_reference) /
                        mmg_N_rdot_reference,
                    ),
                    abs(
                        (exact_prime.N_rdot - mmg_N_rdot_reference) /
                        mmg_N_rdot_reference,
                    ),
                ]
                append!(row, last.(derivative_columns("strict", strict_prime)))
                append!(row, last.(derivative_columns("exact", exact_prime)))
                println(io, join(row, ','))
                write_velocity_comparison_rows(
                    comparison_io,
                    shape,
                    mesh,
                    maximum_section,
                    velocity_mask,
                    (
                        strict_whole_hull=strict_prime,
                        exact_whole_hull=exact_prime,
                        strict_schmitz=strict_schmitz_prime,
                        exact_schmitz=exact_schmitz_prime,
                    ),
                    mmg_velocity,
                )
                println((
                    shape=shape,
                    faces=mesh.nfaces,
                    reciprocity_error=reciprocity_error,
                ))
            end
        end
    end
end

function main()
    default_data_directory = joinpath(@__DIR__, "data")
    data_directory = argument_value("--data-dir=", default_data_directory)
    results_directory = argument_value(
        "--results-dir=",
        joinpath(@__DIR__, "results"),
    )
    geometry_only = "--geometry-only" in ARGS
    fine = "--fine" in ARGS
    mkpath(results_directory)
    source_paths = workshop_paths(data_directory)

    geometry_path = joinpath(results_directory, "kvlcc2_geometry_convergence.csv")
    write_geometry_convergence(geometry_path, source_paths)
    println("Wrote $(geometry_path)")

    if !geometry_only
        derivative_path = joinpath(results_directory, "kvlcc2_derivative_convergence.csv")
        comparison_path = joinpath(
            results_directory,
            "kvlcc2_mmg_velocity_comparison.csv",
        )
        write_derivative_convergence(
            derivative_path,
            comparison_path,
            source_paths;
            fine,
        )
        println("Wrote $(derivative_path)")
        println("Wrote $(comparison_path)")
    end
end

main()
