
function run_dyn_model(
    u::Array{FT, 1},
    u_names::Array{String, 1},
    config::Dict;
    RS::Union{Nothing, ReferenceStatistics} = nothing,
    case_numbers::Vector{Int} = Int[],
)

    if RS === nothing && isempty(case_numbers)
        case_numbers = collect(1:length(config["observations"]["cases"]))
    elseif RS !== nothing && isempty(case_numbers)
        case_numbers = RS.case_numbers
    elseif RS !== nothing && !isempty(case_numbers)
        @warn("Both RS and case_numbers are given! To avoid errors, RS.case_numbers will be used!")
        case_numbers = RS.case_numbers
    end

    if config["model"]["model"] == "KiD"
        sim_vec = run_KiD_multiple_cases(u, u_names, config, case_numbers)
    elseif config["model"]["model"] == "test_model" # for testing
        @assert config["observations"]["data_names"] == ["test data"]
        sim_vec = test_model(u, u_names, config, case_numbers)
    else
        error("Invalid model!")
    end

    outputs = RS === nothing ? sim_vec : outputs = pca_transform(normalize_sim(sim_vec, RS), RS)

    return outputs
end

function run_KiD_multiple_cases(u::Array{FT, 1}, u_names::Array{String, 1}, config::Dict, case_numbers::Vector{Int})

    @assert !isempty(case_numbers)

    outputs = Float64[]
    for case in config["observations"]["cases"][case_numbers]
        config["model"]["w1"] = case.w1
        config["model"]["p0"] = case.p0
        config["model"]["Nd"] = case.Nd
        ode_sol, aux, precip = run_KiD(u, u_names, config["model"])
        single_case_Gvector =
            config["model"]["filter"]["apply"] ?
            ODEsolution2Gvector(ode_sol, aux, precip, config["observations"]["data_names"], config["model"]["filter"]) :
            ODEsolution2Gvector(ode_sol, aux, precip, config["observations"]["data_names"])
        outputs = [outputs; single_case_Gvector]
    end

    return outputs
end

function run_KiD(u::Array{FT, 1}, u_names::Array{String, 1}, model_settings::Dict)

    update_parameters!(model_settings, u, u_names)
    kid_params = create_kid_parameters(FT, model_settings)

    moisture, precip = KD.get_moisture_and_precipitation_types(
        FT,
        model_settings["moisture_choice"],
        model_settings["precipitation_choice"],
        model_settings["rain_formation_choice"],
        model_settings["sedimentation_choice"],
        model_settings["toml_dict"],
    )

    TS = KD.TimeStepping(model_settings["dt"], model_settings["dt_calib"], model_settings["t_end"])
    space, face_space =
        KD.make_function_space(FT, model_settings["z_min"], model_settings["z_max"], model_settings["n_elem"])
    coord = CC.Fields.coordinate_field(space)

    ρ_profile = KD.ρ_ivp(FT, kid_params, model_settings["thermo_params"])
    init = map(coord -> KD.init_1d_column(FT, kid_params, model_settings["thermo_params"], ρ_profile, coord.z), coord)
    Y = KD.initialise_state(moisture, precip, init)
    aux = KD.initialise_aux(
        FT,
        init,
        kid_params,
        model_settings["thermo_params"],
        model_settings["air_params"],
        model_settings["activation_params"],
        TS,
        nothing,
        face_space,
        moisture,
    )
    ode_rhs! = KD.make_rhs_function(moisture, precip)
    problem = ODE.ODEProblem(ode_rhs!, Y, (model_settings["t_ini"], model_settings["t_end"]), aux)
    saveat = model_settings["filter"]["apply"] ? model_settings["filter"]["saveat_t"] : model_settings["t_calib"]
    solution = ODE.solve(problem, ODE.SSPRK33(), dt = model_settings["dt"], saveat = saveat)

    return solution, aux, precip
end

function ODEsolution2Gvector(ODEsol, aux, precip, variables)

    outputs = Float64[]
    for u in ODEsol.u
        for var in variables
            outputs = [outputs; KD.get_variable_data_from_ODE(u, aux, precip, var)]
        end
    end

    return outputs[:]
end

function ODEsolution2Gvector(ODEsol, aux, precip, variables, filter)

    outputs = Float64[]
    _nz_filtered = filter["nz_filtered"]
    _nz_per_filtered_cell = filter["nz_per_filtered_cell"]
    _nz = _nz_filtered * _nz_per_filtered_cell
    _nt_filtered = filter["nt_filtered"]
    _nt_per_filtered_cell = filter["nt_per_filtered_cell"]
    _nvar = length(variables)

    for i in 1:_nt_filtered
        _single_filtered_cell_data = zeros(_nz * _nvar, _nt_per_filtered_cell)
        for j in 1:_nt_per_filtered_cell
            u = ODEsol[(i - 1) * _nt_per_filtered_cell + j]
            for (k, var) in enumerate(variables)
                _single_filtered_cell_data[((k - 1) * _nz + 1):(k * _nz), j] =
                    KD.get_variable_data_from_ODE(u, aux, precip, var)
            end
        end
        single_var_filtered_vec = [
            mean(_single_filtered_cell_data[((j - 1) * _nz_per_filtered_cell + 1):(j * _nz_per_filtered_cell), :])
            for j in 1:(_nz_filtered * _nvar)
        ]
        outputs = [outputs; single_var_filtered_vec]
    end

    return outputs
end

function update_parameters!(model_settings::Dict, u::Array{FT, 1}, u_names::Array{String, 1})
    for (i, name) in enumerate(u_names)
        model_settings["toml_dict"][name]["value"] = u[i]
    end
end

function create_kid_parameters(FT, model_settings::Dict)
    precip_sources = if ("precip_sources" in keys(model_settings))
        model_settings["precip_sources"]
    else
        1
    end
    precip_sinks = if ("precip_sinks" in keys(model_settings))
        model_settings["precip_sinks"]
    else
        1
    end
    r_dry = if ("r_dry" in keys(model_settings))
        model_settings["r_dry"]
    else
        0.04 * 1e-6
    end
    std_dry = if ("std_dry" in keys(model_settings))
        model_settings["std_dry"]
    else
        1.4
    end
    κ = if ("κ" in keys(model_settings))
        model_settings["κ"]
    else
        0.9
    end
    kid_params = KD.Parameters.Kinematic1DParameters{FT}(;
        w1 = model_settings["w1"],
        t1 = model_settings["t1"],
        p0 = model_settings["p0"],
        precip_sources = precip_sources,
        precip_sinks = precip_sinks,
        prescribed_Nd = model_settings["Nd"],
        qtot_flux_correction = Int(model_settings["qtot_flux_correction"]),
        r_dry = r_dry,
        std_dry = std_dry,
        κ = κ,
    )
    return kid_params
end

function test_model(u::Array{FT, 1}, u_names::Array{String, 1}, config::Dict, case_numbers::Vector{Int})
    (n_c, n_v, n_z, n_t) = get_numbers_from_config(config)

    @assert !isempty(case_numbers)
    @assert length(u) == length(u_names) == 2
    @assert Set(u_names) == Set(["a", "b"])
    @assert n_v == 1

    a = u[findall(name -> name == "a", u_names)]
    b = u[findall(name -> name == "b", u_names)]

    n_zt = n_z * n_t
    x = range(0.0, 1.0, n_zt)

    outputs = Float64[]
    for case in config["observations"]["cases"][case_numbers]
        p = case.power
        y = a .* x .^ p .+ b
        outputs = [outputs; y]
    end

    return outputs
end
