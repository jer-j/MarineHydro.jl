using MarineHydro
using TOML

const KVLCC2_REYNOLDS_NUMBER = 4.6e6
const KVLCC2_FROUDE_NUMBER = 0.142
const KVLCC2_DRAFT_RATIO = 20.8 / 320.0
const WATER_DENSITY = 1025.0
const REFERENCE_PATH = joinpath(
    @__DIR__,
    "..",
    "kvlcc2_maneuvering",
    "reference.toml",
)

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

function argument_value(prefix, default)
    argument_index = findfirst(argument -> startswith(argument, prefix), ARGS)
    return isnothing(argument_index) ? default : split(
        ARGS[argument_index],
        "=";
        limit=2,
    )[2]
end

function velocity_reference()
    reference = TOML.parsefile(REFERENCE_PATH)
    values = reference["mmg"]["linear_hull"]
    mmg = MMGLinearHullDerivatives(
        values["Y_v"],
        values["Y_r"],
        values["N_v"],
        values["N_r"],
    )
    return mmg_to_wang_velocity_derivatives(mmg, 20.8, 320.0)
end

velocity_fields(derivatives) = (
    Y_v=derivatives.Y_v,
    Y_r=derivatives.Y_r,
    N_v=derivatives.N_v,
    N_r=derivatives.N_r,
)

function write_method_rows(
    io,
    shape,
    grid,
    method,
    derivatives,
    reference,
    separated_panels,
    linearization,
    coupling_iterations,
    coupling_residual,
    correction_residuals,
    transpiration_flux_derivatives,
)
    for coefficient in keys(reference)
        value = getproperty(derivatives, coefficient)
        target = getproperty(reference, coefficient)
        difference = value - target
        row = (
            "$(shape[1])x$(shape[2])",
            grid.mesh.nfaces,
            length(grid.strips),
            KVLCC2_REYNOLDS_NUMBER,
            method,
            coefficient,
            value,
            target,
            difference,
            abs(difference) / abs(target),
            separated_panels,
            linearization,
            coupling_iterations,
            coupling_residual,
            correction_residuals.sway,
            correction_residuals.yaw,
            transpiration_flux_derivatives[1],
            transpiration_flux_derivatives[2],
        )
        println(io, join(row, ','))
    end
end

function main()
    data_directory = argument_value(
        "--data-dir=",
        joinpath(@__DIR__, "data"),
    )
    output_path = argument_value(
        "--output=",
        joinpath(@__DIR__, "results", "kvlcc2_boundary_layer_comparison.csv"),
    )
    fine = "--fine" in ARGS
    shape_argument = argument_value("--shape=", "")
    shapes = if !isempty(shape_argument)
        dimensions = parse.(Int, split(lowercase(shape_argument), 'x'))
        length(dimensions) == 2 || error("--shape must have the form NIxNJ")
        ((dimensions[1], dimensions[2]),)
    elseif fine
        ((8, 5), (10, 6), (12, 7), (16, 9))
    else
        ((8, 5), (10, 6), (12, 7))
    end
    coupling_iterations = parse(
        Int,
        argument_value("--coupling-iterations=", "0"),
    )
    source_paths = workshop_paths(data_directory)
    reference = velocity_fields(velocity_reference())
    forward_speed = KVLCC2_FROUDE_NUMBER * sqrt(SETTINGS.g)
    kinematic_viscosity = forward_speed / KVLCC2_REYNOLDS_NUMBER
    maximum_section = gothenburg_maximum_section(source_paths)

    mkpath(dirname(output_path))
    open(output_path, "w") do io
        println(
            io,
            join((
                "shape",
                "faces",
                "strips",
                "Re",
                "method",
                "coefficient",
                "value_wang_prime",
                "experimental_mmg_wang_prime",
                "signed_difference",
                "absolute_relative_difference",
                "separated_panels",
                "linearization",
                "coupling_iterations",
                "coupling_residual",
                "transpiration_bem_sway_residual",
                "transpiration_bem_yaw_residual",
                "integrated_dvn_dv",
                "integrated_dvn_dr",
            ), ','),
        )
        for shape in shapes
            grid = read_gothenburg2010_panel_grid(
                source_paths;
                target_shape=shape,
            )
            mesh = grid.mesh
            inviscid = solve_potential_flow_maneuvering(
                mesh,
                forward_speed;
                rho=WATER_DENSITY,
            )
            surge = solve_rigid_body_potential(mesh, :surge)
            correction = viscous_maneuvering_correction(
                grid,
                surge,
                inviscid,
                forward_speed,
                kinematic_viscosity;
                rho=WATER_DENSITY,
                coupling_iterations,
            )
            corrected = apply_viscous_correction(
                inviscid.derivatives,
                correction,
            )
            schmitz_mask = wang_stern_mask(mesh, maximum_section.x)
            schmitz = potential_flow_hydrodynamic_derivatives(
                mesh,
                inviscid.potentials.sway,
                inviscid.potentials.yaw,
                inviscid.potential_gradients.sway[:, 1],
                inviscid.potential_gradients.yaw[:, 1],
                forward_speed;
                rho=WATER_DENSITY,
                velocity_mask=schmitz_mask,
            )

            whole_hull_prime = nondimensionalize_maneuvering_derivatives(
                inviscid.derivatives,
                1.0,
                forward_speed;
                rho=WATER_DENSITY,
            )
            schmitz_prime = nondimensionalize_maneuvering_derivatives(
                schmitz,
                1.0,
                forward_speed;
                rho=WATER_DENSITY,
            )
            pressure_prime = nondimensionalize_viscous_derivatives(
                correction.pressure_derivatives,
                1.0,
                forward_speed;
                rho=WATER_DENSITY,
            )
            shear_prime = nondimensionalize_viscous_derivatives(
                correction.shear_derivatives,
                1.0,
                forward_speed;
                rho=WATER_DENSITY,
            )
            total_correction_prime = nondimensionalize_viscous_derivatives(
                correction.derivatives,
                1.0,
                forward_speed;
                rho=WATER_DENSITY,
            )
            corrected_prime = nondimensionalize_maneuvering_derivatives(
                corrected,
                1.0,
                forward_speed;
                rho=WATER_DENSITY,
            )
            flux_derivatives = vec(sum(
                correction.transpiration_derivatives .* mesh.areas;
                dims=1,
            ))
            separated_panels = count(correction.base_boundary_layer.separated)
            diagnostics = (
                separated_panels,
                correction.linearization,
                correction.coupling_iterations,
                correction.coupling_residual,
                correction.boundary_residuals,
                flux_derivatives,
            )
            methods = (
                whole_hull_inviscid=velocity_fields(whole_hull_prime),
                schmitz_inviscid=velocity_fields(schmitz_prime),
                head_pressure_correction=velocity_fields(pressure_prime),
                head_shear_correction=velocity_fields(shear_prime),
                head_total_correction=velocity_fields(total_correction_prime),
                whole_hull_plus_head=velocity_fields(corrected_prime),
            )
            for (method, values) in pairs(methods)
                write_method_rows(
                    io,
                    shape,
                    grid,
                    method,
                    values,
                    reference,
                    diagnostics...,
                )
            end
            println((
                shape,
                faces=mesh.nfaces,
                separated_panels,
                corrected=methods.whole_hull_plus_head,
            ))
        end
    end
    println("wrote $(output_path)")
end

main()
