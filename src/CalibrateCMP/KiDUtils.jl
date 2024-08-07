
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
    elseif config["model"]["model"] == "KiD_col_sed"
        sim_vec = run_KiD_col_sed_multiple_cases(u, u_names, config, case_numbers)
    elseif config["model"]["model"] == "Box"
        sim_vec = run_box_multiple_cases(u, u_names, config, case_numbers)
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

        if "t_cal" in collect(keys(case))
            config["model"]["filter"] = make_filter_props(
                config["model"]["filter"]["nz_unfiltered"], 
                case.t_cal;
                apply = config["model"]["filter"]["apply"],
                nz_per_filtered_cell = config["model"]["filter"]["nz_per_filtered_cell"],
                nt_per_filtered_cell = config["model"]["filter"]["nt_per_filtered_cell"],
            )
        end
        
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
    common_params = create_common_parameters(
        FT,
        precip_sources = model_settings["precip_sources"],
        precip_sinks = model_settings["precip_sinks"],
        Nd = model_settings["Nd"],
        open_system_activation = model_settings["open_system_activation"],
    )
    kid_params = create_kid_parameters(
        FT,
        w1 = model_settings["w1"],
        t1 = model_settings["t1"],
        p0 = model_settings["p0"],
        qtot_flux_correction = model_settings["qtot_flux_correction"],
        r_dry = model_settings["r_dry"],
        std_dry = model_settings["std_dry"],
        κ = model_settings["κ"],
    )

    moisture = CO.get_moisture_type(model_settings["moisture_choice"], model_settings["toml_dict"])
    precip = CO.get_precipitation_type(
        model_settings["precipitation_choice"],
        model_settings["toml_dict"];
        rain_formation_choice = model_settings["rain_formation_choice"],
        sedimentation_choice = model_settings["sedimentation_choice"],
    )

    TS = CO.TimeStepping(model_settings["dt"], model_settings["dt_calib"], model_settings["t_end"])
    space, face_space =
        KD.make_function_space(FT, model_settings["z_min"], model_settings["z_max"], model_settings["n_elem"])
    coord = CC.Fields.coordinate_field(space)

    ρ_profile = CO.ρ_ivp(FT, kid_params, model_settings["thermo_params"])
    init = map(
        coord -> CO.initial_condition_1d(
            FT,
            common_params,
            kid_params,
            model_settings["thermo_params"],
            ρ_profile,
            coord.z,
        ),
        coord,
    )
    Y = CO.initialise_state(moisture, precip, init)
    aux = KD.initialise_aux(
        FT,
        init,
        common_params,
        kid_params,
        model_settings["thermo_params"],
        model_settings["air_params"],
        model_settings["activation_params"],
        TS,
        nothing,
        face_space,
        moisture,
        precip,
    )
    ode_rhs! = KD.make_rhs_function(moisture, precip)
    problem = ODE.ODEProblem(ode_rhs!, Y, (model_settings["t_ini"], model_settings["t_end"]), aux)
    saveat = model_settings["filter"]["apply"] ? model_settings["filter"]["saveat_t"] : model_settings["t_calib"]
    solution = ODE.solve(problem, ODE.SSPRK33(), dt = model_settings["dt"], saveat = saveat)

    return solution, aux, precip
end

"""
    Run 1D rainshaft simulation with only collisions and sedimentation
"""
function run_KiD_col_sed_multiple_cases(
    u::Array{FT, 1},
    u_names::Array{String, 1},
    config::Dict,
    case_numbers::Vector{Int},
)

    @assert !isempty(case_numbers)

    outputs = Float64[]
    for case in config["observations"]["cases"][case_numbers]
        config["model"]["qt"] = case.qt
        config["model"]["Nd"] = case.Nd
        config["model"]["k"] = case.k
        ode_sol, aux, precip = run_KiD_col_sed(u, u_names, config["model"])
        single_case_Gvector =
            config["model"]["filter"]["apply"] ?
            ODEsolution2Gvector(ode_sol, aux, precip, config["observations"]["data_names"], config["model"]["filter"]) :
            ODEsolution2Gvector(ode_sol, aux, precip, config["observations"]["data_names"])
        outputs = [outputs; single_case_Gvector]
    end

    return outputs
end

function run_KiD_col_sed(u::Array{FT, 1}, u_names::Array{String, 1}, model_settings::Dict)

    update_parameters!(model_settings, u, u_names)
    apply_param_dependency!(model_settings)
    model_settings["toml_dict"]["SB2006_cloud_gamma_distribution_parameter"]["value"] = model_settings["k"]
    common_params = create_common_parameters(
        FT,
        precip_sources = model_settings["precip_sources"],
        precip_sinks = model_settings["precip_sinks"],
        Nd = model_settings["Nd"],
    )
    kid_params = create_kid_parameters(FT)

    moisture = CO.get_moisture_type("NonEquilibriumMoisture", model_settings["toml_dict"])
    precip = CO.get_precipitation_type(
        model_settings["precipitation_choice"],
        model_settings["toml_dict"];
        rain_formation_choice = model_settings["rain_formation_choice"],
        sedimentation_choice = model_settings["sedimentation_choice"],
    )

    TS = CO.TimeStepping(model_settings["dt"], model_settings["dt_calib"], model_settings["t_end"])
    space, face_space =
        KD.make_function_space(FT, model_settings["z_min"], model_settings["z_max"], model_settings["n_elem"])
    coord = CC.Fields.coordinate_field(space)

    init = map(
        coord -> CO.initial_condition_0d(
            FT,
            model_settings["thermo_params"],
            model_settings["qt"],
            model_settings["Nd"],
            model_settings["k"],
            model_settings["rhod"],
        ),
        coord,
    )
    Y = CO.initialise_state(moisture, precip, init)
    aux = KD.initialise_aux(
        FT,
        init,
        common_params,
        kid_params,
        model_settings["thermo_params"],
        model_settings["air_params"],
        model_settings["activation_params"],
        TS,
        nothing,
        face_space,
        moisture,
        precip,
    )
    ode_rhs! = KD.make_rhs_function_col_sed(moisture, precip)
    problem = ODE.ODEProblem(ode_rhs!, Y, (model_settings["t_ini"], model_settings["t_end"]), aux)
    saveat = model_settings["filter"]["apply"] ? model_settings["filter"]["saveat_t"] : model_settings["t_calib"]
    solution = ODE.solve(problem, ODE.SSPRK33(), dt = model_settings["dt"], saveat = saveat)

    return solution, aux, precip
end

"""
    Run 0D box simulations with only collisional processes
"""
function run_box_multiple_cases(u::Array{FT, 1}, u_names::Array{String, 1}, config::Dict, case_numbers::Vector{Int})

    @assert !isempty(case_numbers)

    outputs = Float64[]
    for case in config["observations"]["cases"][case_numbers]
        config["model"]["qt"] = case.qt
        config["model"]["Nd"] = case.Nd
        config["model"]["k"] = case.k
        ode_sol, aux, precip = run_box(u, u_names, config["model"])
        single_case_Gvector =
            config["model"]["filter"]["apply"] ?
            ODEsolution2Gvector(ode_sol, aux, precip, config["observations"]["data_names"], config["model"]["filter"]) :
            ODEsolution2Gvector(ode_sol, aux, precip, config["observations"]["data_names"])
        outputs = [outputs; single_case_Gvector]
    end

    return outputs
end

function run_box(u::Array{FT, 1}, u_names::Array{String, 1}, model_settings::Dict)

    update_parameters!(model_settings, u, u_names)
    apply_param_dependency!(model_settings)
    model_settings["toml_dict"]["SB2006_cloud_gamma_distribution_parameter"]["value"] = model_settings["k"]
    common_params = create_common_parameters(
        FT,
        precip_sources = model_settings["precip_sources"],
        precip_sinks = model_settings["precip_sinks"],
        Nd = model_settings["Nd"],
    )

    moisture = CO.get_moisture_type("NonEquilibriumMoisture", model_settings["toml_dict"])
    precip = CO.get_precipitation_type(
        model_settings["precipitation_choice"],
        model_settings["toml_dict"];
        rain_formation_choice = model_settings["rain_formation_choice"],
    )

    TS = CO.TimeStepping(model_settings["dt"], model_settings["dt_calib"], model_settings["t_end"])

    domain = CC.Domains.IntervalDomain(
        CC.Geometry.ZPoint{FT}(0),
        CC.Geometry.ZPoint{FT}(1),
        boundary_names = (:bottom, :top),
    )
    mesh = CC.Meshes.IntervalMesh(domain, nelems = 1)
    space = CC.Spaces.CenterFiniteDifferenceSpace(mesh)
    coord = CC.Fields.coordinate_field(space)

    init = map(
        coord -> CO.initial_condition_0d(
            FT,
            model_settings["thermo_params"],
            model_settings["qt"],
            model_settings["Nd"],
            model_settings["k"],
            model_settings["rhod"],
        ),
        coord,
    )
    Y = CO.initialise_state(moisture, precip, init)
    aux = CO.initialise_aux(
        FT,
        init,
        common_params,
        model_settings["thermo_params"],
        model_settings["air_params"],
        model_settings["activation_params"],
        TS,
        nothing,
        moisture,
        precip,
    )
    ode_rhs! = BX.make_rhs_function(moisture, precip)
    problem = ODE.ODEProblem(ode_rhs!, Y, (model_settings["t_ini"], model_settings["t_end"]), aux)
    saveat = model_settings["filter"]["apply"] ? model_settings["filter"]["saveat_t"] : model_settings["t_calib"]
    solution = ODE.solve(problem, ODE.SSPRK33(), dt = model_settings["dt"], saveat = saveat)

    return solution, aux, precip
end

"""
    General functions for converting ode results to vector of values
"""
function ODEsolution2Gvector(ODEsol, aux, precip, variables)

    outputs = Float64[]
    for u in ODEsol.u
        for var in variables
            outputs = [outputs; CO.get_variable_data_from_ODE(u, aux, precip, var)]
        end
    end

    return outputs[:]
end

function ODEsolution2Gvector(ODEsol, aux, precip, variables, filter)

    outputs = Float64[]
    _nz_filtered = filter["nz_filtered"]
    _nz_per_filtered_cell = filter["nz_per_filtered_cell"]
    _nz = _nz_filtered .* _nz_per_filtered_cell
    _nt_filtered = filter["nt_filtered"]
    _nt_per_filtered_cell = filter["nt_per_filtered_cell"]
    _n_variables = length(variables)

    _n_single_time_unfiltered = sum(_nz)

    for i in 1:_nt_filtered
        # get unfiltered data in a filtered time interval 
        _single_filtered_cell_data = zeros(_n_single_time_unfiltered, _nt_per_filtered_cell)
        for j in 1:_nt_per_filtered_cell
            u = ODEsol[:, (i - 1) * _nt_per_filtered_cell + j]
            for (k, var) in enumerate(variables)
                _range = (sum(_nz[1:(k - 1)]) + 1):sum(_nz[1:k])
                _single_filtered_cell_data[_range, j] = CO.get_variable_data_from_ODE(u, aux, precip, var)
            end
        end

        # filter data spatially
        for j in 1:_n_variables
            for k in 1:_nz_filtered[j]
                _ind = sum(_nz[1:(j - 1)]) + (k - 1) * _nz_per_filtered_cell[j]
                _range = (_ind + 1):(_ind + _nz_per_filtered_cell[j])
                outputs = [outputs; mean(_single_filtered_cell_data[_range, :])]
            end
        end

    end

    return outputs
end

function update_parameters!(model_settings::Dict, u::Array{FT, 1}, u_names::Array{String, 1})
    for (i, name) in enumerate(u_names)
        model_settings["toml_dict"][name]["value"] = u[i]
    end
end

function apply_param_dependency!(model_settings::Dict)
    if "param_dependencies" in keys(model_settings)
        for dep in model_settings["param_dependencies"]
            model_settings["toml_dict"][dep.dependant]["value"] =
                model_settings["toml_dict"][dep.base]["value"] * dep.ratio
        end
    end
end

function create_common_parameters(FT; precip_sources = 1, precip_sinks = 1, Nd = 1e8, open_system_activation = false)
    common_params = CO.Parameters.CommonParameters{FT}(;
        precip_sources = precip_sources,
        precip_sinks = precip_sinks,
        prescribed_Nd = Nd,
        open_system_activation = open_system_activation,
    )
    return common_params
end

function create_kid_parameters(
    FT;
    w1 = 2.0,
    t1 = 600,
    p0 = 100000,
    z_0 = 0.0,
    z_1 = 740.0,
    z_2 = 3260.0,
    rv_0 = 0.015,
    rv_1 = 0.0138,
    rv_2 = 0.0024,
    tht_0 = 297.9,
    tht_1 = 297.9,
    tht_2 = 312.66,
    qtot_flux_correction = 1,
    r_dry = 0.04 * 1e-6,
    std_dry = 1.4,
    κ = 0.9,
)
    kid_params = KD.Parameters.KinematicDriverParameters{FT}(;
        w1,
        t1,
        p0,
        z_0,
        z_1,
        z_2,
        rv_0,
        rv_1,
        rv_2,
        tht_0,
        tht_1,
        tht_2,
        qtot_flux_correction,
        r_dry,
        std_dry,
        κ,
    )
    return kid_params
end

function test_model(u::Array{FT, 1}, u_names::Array{String, 1}, config::Dict, case_numbers::Vector{Int})
    (; n_heights, n_times) = get_numbers_from_config(config)

    @assert !isempty(case_numbers)
    @assert length(u) == length(u_names) == 2
    @assert Set(u_names) == Set(["a", "b"])

    a = u[findall(name -> name == "a", u_names)]
    b = u[findall(name -> name == "b", u_names)]

    n_single_case = sum(n_heights) * n_times
    x = range(0.0, 1.0, n_single_case)

    outputs = Float64[]
    for case in config["observations"]["cases"][case_numbers]
        p = case.power
        y = a .* x .^ p .+ b
        outputs = [outputs; y]
    end

    return outputs
end
