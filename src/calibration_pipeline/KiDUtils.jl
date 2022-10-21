
function run_dyn_model(
    u::Array{FT, 1},
    u_names::Array{String, 1},
    config::Dict;
    RS::Union{Nothing, ReferenceStatistics} = nothing,
)

    if config["model"]["model"] == "KiD"
        outputs = run_KiD_multiple_cases(u, u_names, config)
    elseif config["model"]["model"] == "terminal_velocity"
        @assert config["observations"]["data_names"] == ["rho", "qr", "rain averaged terminal velocity"]
        @assert RS != nothing
        outputs = compute_terminal_velocity(u, u_names, config, RS)
    else
        error("Invalid model!")
    end

    if RS != nothing
        outputs = RS.P_pca' * outputs
    end

    return outputs
end

function run_KiD_multiple_cases(u::Array{FT, 1}, u_names::Array{String, 1}, config::Dict)

    outputs = Float64[]
    for (i, case) in enumerate(config["observations"]["cases"])
        config["model"]["w1"] = case.w1
        config["model"]["p0"] = case.p0
        config["model"]["Nd"] = case.Nd
        ode_sol, aux = run_KiD(u, u_names, config["model"])
        outputs = [
            outputs
            ODEsolution2Gvector(
                ode_sol,
                aux,
                config["observations"]["data_names"],
                config["observations"]["ynorm"][i, :],
            )
        ]
    end

    return outputs
end

function run_KiD(u::Array{FT, 1}, u_names::Array{String, 1}, model_settings::Dict)

    params_calib = create_param_dict(u, u_names)
    params = create_parameter_set(FT, model_settings, params_calib)

    moisture, precip = equation_types(
        model_settings["moisture_choice"],
        model_settings["precipitation_choice"],
        model_settings["rain_formation_choice"],
    )

    TS = TimeStepping(model_settings["dt"], model_settings["dt_output"], model_settings["t_end"])
    space, face_space =
        make_function_space(FT, model_settings["z_min"], model_settings["z_max"], model_settings["n_elem"])
    coord = CC.Fields.coordinate_field(space)

    ρ_profile = ρ_ivp(FT, params)
    init = map(coord -> init_1d_column(FT, params, ρ_profile, coord.z), coord)
    Y = initialise_state(moisture, precip, init)
    aux = initialise_aux(FT, init, params, TS, nothing, face_space, moisture)
    ode_rhs! = make_rhs_function(moisture, precip)
    problem = ODE.ODEProblem(ode_rhs!, Y, (model_settings["t_ini"], model_settings["t_end"]), aux)
    solution = ODE.solve(problem, ODE.SSPRK33(), dt = model_settings["dt"], saveat = model_settings["t_calib"])

    return solution, aux
end

function ODEsolution2Gvector(ODEsol, aux, variables, norm_vec)

    outputs = Float64[]

    ρ_dry = parent(aux.moisture_variables.ρ_dry)
    for u in ODEsol.u
        ρ = ρ_dry .+ parent(u.ρq_tot)

        for (i, var) in enumerate(variables)
            outputs = [outputs; get_variable_data_from_ODE(u, ρ[:], var) ./ norm_vec[i]]
        end
    end

    return outputs
end

function get_variable_data_from_ODE(u, ρ::Vector{Float64}, var::String)

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
    else
        error("Data name \"" * var * "\" not recognized!!")
    end

    return output

end

function compute_terminal_velocity(u::Array{FT, 1}, u_names::Array{String, 1}, config::Dict, RS::ReferenceStatistics)
    return compute_terminal_velocity(u, u_names, config, RS.y_full)
end

Base.broadcastable(ps::CM.CommonTypes.RainType) = Ref(ps)
Base.broadcastable(ps::CM.Parameters.CloudMicrophysicsParameters) = Ref(ps)
function compute_terminal_velocity(u::Array{FT, 1}, u_names::Array{String, 1}, config::Dict, obs_vec_full::Array{FT, 1})

    model_settings = config["model"]
    params_calib = create_param_dict(u, u_names)
    params = create_parameter_set(FT, model_settings, params_calib)
    microphys_params = Parameters.microphysics_params(params)

    _n_elem = model_settings["n_elem"]
    _n_times = length(model_settings["t_calib"])
    _n_var = 3
    _n_single_sim = _n_elem * _n_times * _n_var
    _n_cases = Int(length(obs_vec_full) / _n_single_sim)

    outputs = Float64[]
    for i in 1:_n_cases
        y_vec = obs_vec_full[((i - 1) * _n_single_sim + 1):(i * _n_single_sim)]
        y_mat = reshape(y_vec, :, _n_times)
        ρ_mat = y_mat[1:_n_elem, :]
        q_rai_mat = y_mat[(_n_elem + 1):(2 * _n_elem), :]

        norm_vec = config["observations"]["ynorm"][i, :]

        for j in 1:_n_times
            ρ = ρ_mat[:, j] .* norm_vec[1]
            q_rai = q_rai_mat[:, j] .* norm_vec[2]
            term_vel_rai = CM1.terminal_velocity.(microphys_params, CM.CommonTypes.RainType(), ρ, q_rai)
            outputs = [outputs; ρ ./ norm_vec[1]; q_rai ./ norm_vec[2]; term_vel_rai ./ norm_vec[3]]
        end

    end
    return outputs
end

# Equations to solve for mositure and precipitation variables
function equation_types(moisture_choice::String, precipitation_choice::String, rain_formation_choice::String)
    if moisture_choice == "EquilibriumMoisture"
        moisture = EquilibriumMoisture_ρdTq()
    elseif moisture_choice == "NonEquilibriumMoisture"
        moisture = NonEquilibriumMoisture_ρdTq()
    else
        error("Invalid moisture choice: $moisture_choice")
    end
    if precipitation_choice == "NoPrecipitation"
        precip = NoPrecipitation()
    elseif precipitation_choice == "Precipitation0M"
        precip = Precipitation0M()
    elseif precipitation_choice == "Precipitation1M"
        if rain_formation_choice == "CliMA_1M"
            precip = Precipitation1M(OneMomentRainFormation())
        elseif rain_formation_choice == "KK2000"
            precip = Precipitation1M(CMT.KK2000Type())
        elseif rain_formation_choice == "B1994"
            precip = Precipitation1M(CMT.B1994Type())
        elseif rain_formation_choice == "TC1980"
            precip = Precipitation1M(CMT.TC1980Type())
        elseif rain_formation_choice == "LD2004"
            precip = Precipitation1M(CMT.LD2004Type())
        else
            error("Invalid rain formation choice: $rain_formation_choice")
        end
    else
        error("Invalid precipitation choice: $precipitation_choice")
    end
    return moisture, precip
end

function create_param_dict(u::Array{FT, 1}, u_names::Array{String, 1}) where {FT <: AbstractFloat}
    params = Dict{String, Any}()
    for (i, name) in enumerate(u_names)
        params[name] = u[i]
    end
    return params
end

function create_parameter_set(FT, model_settings::Dict, params_cal::Dict)
    thermo_params = model_settings["thermo_params"]
    TP = typeof(thermo_params)

    fixed_microphys_pairs = model_settings["fixed_microphys_param_pairs"]
    pairs = []
    for (key, val) in params_cal
        push!(pairs, Pair(Symbol(key), val))
    end
    microphys_params =
        CM.Parameters.CloudMicrophysicsParameters{FT, TP}(; fixed_microphys_pairs..., pairs..., thermo_params)
    MP = typeof(microphys_params)

    param_set = KP.KinematicParameters{FT, MP}(;
        w1 = model_settings["w1"],
        t1 = model_settings["t1"],
        p0 = model_settings["p0"],
        precip_sources = 1,
        precip_sinks = 1,
        prescribed_Nd = model_settings["Nd"],
        qtot_flux_correction = Int(model_settings["qtot_flux_correction"]),
        microphys_params,
    )
    return param_set
end
