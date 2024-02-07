
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

    moisture, precip = equation_types(
        model_settings["moisture_choice"],
        model_settings["precipitation_choice"],
        model_settings["rain_formation_choice"],
        model_settings["sedimentation_choice"],
        model_settings["toml_dict"],
    )

    TS = TimeStepping(model_settings["dt"], model_settings["dt_calib"], model_settings["t_end"])
    space, face_space =
        make_function_space(FT, model_settings["z_min"], model_settings["z_max"], model_settings["n_elem"])
    coord = CC.Fields.coordinate_field(space)

    ρ_profile = ρ_ivp(FT, kid_params, model_settings["thermo_params"])
    init = map(coord -> init_1d_column(FT, kid_params, model_settings["thermo_params"], ρ_profile, coord.z), coord)
    Y = initialise_state(moisture, precip, init)
    aux = initialise_aux(
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
    ode_rhs! = make_rhs_function(moisture, precip)
    problem = ODE.ODEProblem(ode_rhs!, Y, (model_settings["t_ini"], model_settings["t_end"]), aux)
    saveat = model_settings["filter"]["apply"] ? model_settings["filter"]["saveat_t"] : model_settings["t_calib"]
    solution = ODE.solve(problem, ODE.SSPRK33(), dt = model_settings["dt"], saveat = saveat)

    return solution, aux, precip
end

function ODEsolution2Gvector(ODEsol, aux, precip, variables)

    outputs = Float64[]
    for u in ODEsol.u
        for var in variables
            outputs = [outputs; get_variable_data_from_ODE(u, aux, precip, var)]
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
                    get_variable_data_from_ODE(u, aux, precip, var)
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

function get_variable_data_from_ODE(u, aux, precip, var::String)

    ρ = parent(aux.moisture_variables.ρ_dry)[:] .+ parent(u.ρq_tot)[:]

    if var == "qt"
        output = parent(u.ρq_tot) ./ ρ
    elseif var == "ql"
        output = parent(u.ρq_liq) ./ ρ
    elseif var == "qr"
        output = parent(u.ρq_rai) ./ ρ
    elseif var == "qv"
        output = (parent(u.ρq_tot) .- parent(u.ρq_liq)) ./ ρ
    elseif var == "qlr"
        output = (parent(u.ρq_liq) .+ parent(u.ρq_rai)) ./ ρ
    elseif var == "rt"
        _qtot = parent(u.ρq_tot) ./ ρ
        output = _qtot ./ (1 .- _qtot)
    elseif var == "rl"
        _qtot = parent(u.ρq_tot) ./ ρ
        output = parent(u.ρq_liq) ./ ρ ./ (1 .- _qtot)
    elseif var == "rr"
        _qtot = parent(u.ρq_tot) ./ ρ
        output = parent(u.ρq_rai) ./ ρ ./ (1 .- _qtot)
    elseif var == "rv"
        _qtot = parent(u.ρq_tot) ./ ρ
        output = (parent(u.ρq_tot) .- parent(u.ρq_liq)) ./ ρ ./ (1 .- _qtot)
    elseif var == "rlr"
        _qtot = parent(u.ρq_tot) ./ ρ
        output = (parent(u.ρq_liq) .+ parent(u.ρq_rai)) ./ ρ ./ (1 .- _qtot)
    elseif var == "Nl"
        output = parent(u.N_liq)
    elseif var == "Nr"
        output = parent(u.N_rai)
    elseif var == "Na"
        output = parent(u.N_aer)
    elseif var == "rho"
        output = ρ
    elseif var == "rain averaged terminal velocity"
        qr = parent(u.ρq_rai) ./ ρ
        if precip isa Precipitation1M
            vt = CM1.terminal_velocity.(precip.rain, precip.sedimentation.rain, ρ, qr)
        elseif precip isa Precipitation2M
            Nr = parent(u.N_rai)
            ff(xx, yy, zz) = CM2.rain_terminal_velocity(precip.rain_formation, precip.sedimentation, xx, yy, zz)[2]
            vt = ff.(qr, ρ, Nr)
        else
            error("Computing rain averaged terminal velocity for the given precipitation style is invalid!!")
        end
        output = vt
    elseif var == "rainrate"
        qr = parent(u.ρq_rai) ./ ρ
        if precip isa Precipitation1M
            vt = CM1.terminal_velocity.(precip.rain, precip.sedimentation.rain, ρ, qr)
        elseif precip isa Precipitation2M
            Nr = parent(u.N_rai)
            ff(xx, yy, zz) = CM2.rain_terminal_velocity(precip.rain_formation, precip.sedimentation, xx, yy, zz)[2]
            vt = ff.(qr, ρ, Nr)
        else
            error("Computing rainrate for the given precipitation style is invalid!!")
        end
        output = qr .* ρ .* vt .* 3600
    else
        error("Data name \"" * var * "\" not recognized!!")
    end

    return output

end

# Equations to solve for mositure and precipitation variables
function equation_types(
    moisture_choice::String,
    precipitation_choice::String,
    rain_formation_choice::String,
    sedimentation_choice::String,
    toml_dict,
)
    if moisture_choice == "EquilibriumMoisture"
        moisture = EquilibriumMoisture_ρdTq()
    elseif moisture_choice == "NonEquilibriumMoisture"
        moisture = NonEquilibriumMoisture_ρdTq(CMP.CloudLiquid(FT, toml_dict), CMP.CloudIce(FT, toml_dict))
    else
        error("Invalid moisture choice: $moisture_choice")
    end
    if precipitation_choice == "NoPrecipitation"
        precip = NoPrecipitation()
    elseif precipitation_choice == "Precipitation0M"
        precip = Precipitation0M(CMP.Parameters0M(FT, toml_dict))
    elseif precipitation_choice == "Precipitation1M"
        if sedimentation_choice == "CliMA_1M"
            st = CMP.Blk1MVelType(FT, toml_dict)
        elseif sedimentation_choice == "Chen2022"
            st = CMP.Chen2022VelType(FT, toml_dict)
        else
            error("Invalid sedimentation choice: $sedimentation_choice")
        end
        if rain_formation_choice == "CliMA_1M"
            rain_params = CMP.Rain(FT, toml_dict)
            rf = rain_params.acnv1M
        elseif rain_formation_choice == "KK2000"
            rf = CMP.KK2000(FT, toml_dict)
        elseif rain_formation_choice == "B1994"
            rf = CMP.B1994(FT, toml_dict)
        elseif rain_formation_choice == "TC1980"
            rf = CMP.TC1980(FT, toml_dict)
        elseif rain_formation_choice == "LD2004"
            rf = CMP.LD2004(FT, toml_dict)
        elseif rain_formation_choice == "VarTimeScaleAcnv"
            rf = CMP.VarTimescaleAcnv(FT, toml_dict)
        else
            error("Invalid rain formation choice: $rain_formation_choice")
        end
        precip = Precipitation1M(
            CMP.CloudLiquid(FT, toml_dict),
            CMP.CloudIce(FT, toml_dict),
            CMP.Rain(FT, toml_dict),
            CMP.Snow(FT, toml_dict),
            CMP.CollisionEff(FT, toml_dict),
            rf,
            st,
        )
    elseif precipitation_choice == "Precipitation2M"
        if sedimentation_choice == "Chen2022"
            st = CMP.Chen2022VelTypeRain(FT, toml_dict)
        elseif sedimentation_choice == "SB2006"
            st = CMP.SB2006VelType(FT, toml_dict)
        else
            error("Invalid sedimentation choice: $sedimentation_choice")
        end
        if rain_formation_choice == "SB2006"
            precip = Precipitation2M(CMP.SB2006(FT, toml_dict), st)
        else
            error("Invalid rain formation choice: $rain_formation_choice")
        end
    else
        error("Invalid precipitation choice: $precipitation_choice")
    end
    return moisture, precip
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
    kid_params = KP.KinematicParameters{FT}(;
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
