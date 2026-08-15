using MarineHydro
using TOML

const REFERENCE_PATH = joinpath(@__DIR__, "references.toml")

"""
    argument_value(prefix, default)

Return the command-line value following `prefix`, or `default` when absent.
"""
function argument_value(prefix, default)
    index = findfirst(argument -> startswith(argument, prefix), ARGS)
    return isnothing(index) ? default : split(ARGS[index], "="; limit=2)[2]
end

function reference_value(reference, name)
    value = get(reference, name, "NaN")
    return value isa Real ? Float64(value) : NaN
end

function require_file(path)
    isfile(path) || error(
        "Missing public hull geometry $(path). Run " *
        "validation/public_hulls/fetch_geometry.sh first.",
    )
    return path
end

function gothenburg_paths(gothenburg_data, hull)
    root = joinpath(gothenburg_data, hull)
    filenames = hull == "KCS" ?
        ("kcs_bow1.dat", "kcs_stn1.dat") :
        ("kvlcc_bow1.dat", "kvlcc2_stn1.dat")
    return collect(require_file.(joinpath.(root, filenames)))
end

"""
    public_hull_cases(; fine=false, data_directory, gothenburg_data)

Construct KCS, KVLCC2, DTMB 5415, DTC, and Wigley meshes in MarineHydro's
``x``-forward, ``y``-port, ``z``-up convention. `fine=false` keeps the dense
boundary-element study below 500 panels per hull. The fine discretization is
intended for notebook figures and convergence checks.
"""
function public_hull_cases(;
    fine::Bool=false,
    data_directory=joinpath(@__DIR__, "data"),
    gothenburg_data=joinpath(@__DIR__, "..", "gothenburg2010", "data"),
)
    references = TOML.parsefile(REFERENCE_PATH)
    gothenburg_shape = fine ? (18, 10) : (14, 8)
    dtmb_shape = fine ? (37, 13) : (29, 11)
    dtc_clustering = fine ? (24, 8, 5) : (18, 6, 4)
    wigley_shape = fine ? (32, 8) : (24, 6)

    kcs_reference = references["KCS"]
    kcs_mesh = read_gothenburg2010_mesh(
        gothenburg_paths(gothenburg_data, "KCS");
        length_scale=kcs_reference["length_m"],
        target_shape=gothenburg_shape,
    )

    kvlcc2_reference = references["KVLCC2"]
    kvlcc2_mesh = read_gothenburg2010_mesh(
        gothenburg_paths(gothenburg_data, "KVLCC2");
        length_scale=kvlcc2_reference["length_m"],
        target_shape=gothenburg_shape,
    )

    dtmb_reference = references["DTMB5415"]
    dtmb_surface = only(read_plot3d_structured_surfaces(require_file(joinpath(
        data_directory,
        "DTMB5415",
        "5415_static.net",
    ))))
    dtmb_length = dtmb_reference["length_m"]
    dtmb_draft = dtmb_reference["draft_m"]
    dtmb_mesh = structured_surface_waterline_mesh(
        dtmb_surface;
        target_shape=dtmb_shape,
        mirror=true,
        coordinate_transform=point -> [
            point[1] - dtmb_length / 2,
            -point[2],
            point[3] - dtmb_draft,
        ],
        interior_point=[0.0, 0.0, -dtmb_draft / 2],
    )

    dtc_reference = references["DTC"]
    dtc_length = dtc_reference["length_m"]
    dtc_draft = dtc_reference["draft_m"]
    dtc_mesh = read_ascii_stl_waterline_mesh(
        require_file(joinpath(data_directory, "DTC", "DTC-scaled.stl"));
        clustering=dtc_clustering,
        coordinate_transform=point -> [
            point[1] - dtc_length / 2,
            -point[2],
            point[3] - dtc_draft,
        ],
        interior_point=[0.0, 0.0, -dtc_draft / 2],
    )

    wigley_reference = references["Wigley"]
    wigley_mesh = wigley_hull_mesh(
        wigley_reference["length_m"],
        wigley_reference["beam_m"],
        wigley_reference["draft_m"];
        longitudinal_panels=wigley_shape[1],
        vertical_panels=wigley_shape[2],
    )

    meshes = Dict(
        "KCS" => kcs_mesh,
        "KVLCC2" => kvlcc2_mesh,
        "DTMB5415" => dtmb_mesh,
        "DTC" => dtc_mesh,
        "Wigley" => wigley_mesh,
    )
    labels = Dict(
        "KCS" => "KCS",
        "KVLCC2" => "KVLCC2",
        "DTMB5415" => "DTMB 5415",
        "DTC" => "DTC",
        "Wigley" => "Wigley",
    )
    return [
        (
            id,
            name=labels[id],
            mesh=meshes[id],
            reference=references[id],
        )
        for id in ("KCS", "KVLCC2", "DTMB5415", "DTC", "Wigley")
    ]
end

function relative_difference(value, reference)
    return isfinite(reference) && !iszero(reference) ?
        (value - reference) / reference : NaN
end

"""
    analyze_public_hull(case; density=1025, forward_speed=1)

Solve independent sway and yaw potentials, apply the configured Schmitz stern
cutoff only to velocity-dependent integrals, and return geometry, derivative,
and added-mass comparison metrics. Acceleration derivatives always use the
complete wetted surface.
"""
function analyze_public_hull(
    case;
    density::Real=1025.0,
    forward_speed::Real=1.0,
)
    reference = case.reference
    length = reference["length_m"]
    draft = reference["draft_m"]
    displacement = reference["displacement_m3"]
    x_cut = reference["schmitz_x_over_length"] * length
    velocity_mask = wang_stern_mask(case.mesh, x_cut)
    exact = solve_potential_flow_maneuvering(
        case.mesh,
        forward_speed;
        rho=density,
        velocity_mask,
    )
    strict_derivatives = wang_hydrodynamic_derivatives(
        case.mesh,
        exact.potentials.sway,
        exact.potential_gradients.sway[:, 1],
        forward_speed;
        rho=density,
        velocity_mask,
    )
    exact_prime = nondimensionalize_maneuvering_derivatives(
        exact.derivatives,
        length,
        forward_speed;
        rho=density,
    )
    strict_prime = nondimensionalize_maneuvering_derivatives(
        strict_derivatives,
        length,
        forward_speed;
        rho=density,
    )

    added_sway = -exact.derivatives.Y_vdot
    added_yaw = -exact.derivatives.N_rdot
    strict_added_yaw = -strict_derivatives.N_rdot
    A22_over_rho_displacement = added_sway / (density * displacement)
    A66_over_rho_displacement_length2 =
        added_yaw / (density * displacement * length^2)
    m_y_prime = added_sway / (0.5 * density * length^2 * draft)
    J_z_prime = added_yaw / (0.5 * density * length^4 * draft)
    strict_J_z_prime = strict_added_yaw /
        (0.5 * density * length^4 * draft)
    Y_beta_mmg = -exact_prime.Y_v * length / draft
    N_beta_mmg = -exact_prime.N_v * length / draft
    strict_Y_beta_mmg = -strict_prime.Y_v * length / draft
    strict_N_beta_mmg = -strict_prime.N_v * length / draft
    reference_Y_beta = reference_value(reference, "reference_Y_beta_mmg")
    reference_N_beta = reference_value(reference, "reference_N_beta_mmg")
    reference_m_y = reference_value(reference, "reference_m_y_prime")
    reference_J_z = reference_value(reference, "reference_J_z_prime")
    reference_A22 = reference_value(
        reference,
        "reference_A22_over_rho_displacement",
    )
    reference_area = reference["wetted_area_m2"]
    computed_volume = mesh_signed_volume(case.mesh)
    computed_area = mesh_surface_area(case.mesh)
    reciprocity_scale = max(
        abs(exact.derivatives.Y_rdot),
        abs(exact.derivatives.N_vdot),
        eps(Float64),
    )

    return (
        hull=case.name,
        panels=case.mesh.nfaces,
        velocity_panels=count(velocity_mask),
        x_cut_over_length=x_cut / length,
        surface_area=computed_area,
        surface_area_relative_error=reference_area > 0 ?
            computed_area / reference_area - 1 : NaN,
        volume=computed_volume,
        volume_relative_error=computed_volume / displacement - 1,
        sway_boundary_residual=exact.boundary_residuals.sway,
        yaw_boundary_residual=exact.boundary_residuals.yaw,
        acceleration_reciprocity_error=abs(
            exact.derivatives.Y_rdot - exact.derivatives.N_vdot,
        ) / reciprocity_scale,
        A22_over_rho_displacement,
        reference_A22_over_rho_displacement=reference_A22,
        A22_relative_difference=relative_difference(
            A22_over_rho_displacement,
            reference_A22,
        ),
        A66_over_rho_displacement_length2,
        m_y_prime,
        reference_m_y_prime=reference_m_y,
        m_y_relative_difference=relative_difference(m_y_prime, reference_m_y),
        J_z_prime,
        strict_J_z_prime,
        reference_J_z_prime=reference_J_z,
        J_z_relative_difference=relative_difference(J_z_prime, reference_J_z),
        strict_J_z_relative_difference=(strict_J_z_prime - J_z_prime) / J_z_prime,
        Y_beta_mmg,
        strict_Y_beta_mmg,
        reference_Y_beta_mmg=reference_Y_beta,
        Y_beta_relative_difference=relative_difference(Y_beta_mmg, reference_Y_beta),
        N_beta_mmg,
        strict_N_beta_mmg,
        reference_N_beta_mmg=reference_N_beta,
        N_beta_relative_difference=relative_difference(N_beta_mmg, reference_N_beta),
    )
end

function write_results(path, rows)
    isempty(rows) && throw(ArgumentError("at least one result row is required"))
    mkpath(dirname(path))
    columns = propertynames(first(rows))
    open(path, "w") do io
        println(io, join(columns, ','))
        for row in rows
            println(io, join((getproperty(row, column) for column in columns), ','))
        end
    end
    return path
end

function main()
    fine = "--fine" in ARGS
    data_directory = argument_value(
        "--data-dir=",
        joinpath(@__DIR__, "data"),
    )
    gothenburg_data = argument_value(
        "--gothenburg-data=",
        joinpath(@__DIR__, "..", "gothenburg2010", "data"),
    )
    output_path = argument_value(
        "--output=",
        joinpath(
            @__DIR__,
            "results",
            fine ? "public_hull_comparison_fine.csv" : "public_hull_comparison.csv",
        ),
    )
    cases = public_hull_cases(; fine, data_directory, gothenburg_data)
    rows = NamedTuple[]
    for case in cases
        println((hull=case.name, panels=case.mesh.nfaces, status="solving"))
        row = analyze_public_hull(case)
        push!(rows, row)
        println(row)
    end
    write_results(output_path, rows)
    println((output_path, hulls=length(rows)))
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
