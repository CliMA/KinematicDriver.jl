import OrdinaryDiffEq as ODE
import ClimaCore as CC
import ClimaParams as CP
import CloudMicrophysics.Parameters as CMP
import CloudMicrophysics as CM
import KinematicDriver
import KinematicDriver.Common as CO
import KinematicDriver.K1DModel as K1D

include(joinpath(pkgdir(KinematicDriver), "test", "create_parameters.jl"))
include(joinpath(pkgdir(KinematicDriver), "test", "plotting_utils.jl"))

function run_KiD_simulation(::Type{FT}, opts) where {FT}
    if opts["z_2"] < opts["z_max"]
        throw(ArgumentError("z_max cannot be higher than z_2"))
    end

    # Equations to solve for mositure and precipitation variables
    moisture_choice = opts["moisture_choice"]
    precipitation_choice = opts["precipitation_choice"]
    rain_formation_choice = opts["rain_formation_scheme_choice"]
    sedimentation_choice = opts["sedimentation_scheme_choice"]
    @info moisture_choice, precipitation_choice, rain_formation_choice, sedimentation_choice

    # Decide the output folder name based on options
    folder = "Output_$(moisture_choice)_$(precipitation_choice)"
    if precipitation_choice in ["Precipitation1M", "Precipitation2M"]
        folder *= "_$(rain_formation_choice)"
        sedimentation_choice == "Chen2022" && (folder *= "_Chen2022")
    elseif precipitation_choice == "CloudyPrecip"
        folder *= "_$(opts["num_moments"])"
    end
    opts["qtot_flux_correction"] && (folder *= "_wFC")
    opts["open_system_activation"] && (folder *= "_OSA")
    opts["local_activation"] && (folder *= "_LA")
    path = joinpath(@__DIR__, folder)
    mkpath(path)

    # Overwrite the defaults parameters based on options
    default_toml_dict = CP.create_toml_dict(FT)
    toml_dict = override_toml_dict(
        default_toml_dict;
        w1 = FT(opts["w1"]),
        t1 = FT(opts["t1"]),
        p0 = FT(opts["p0"]),
        z_0 = FT(opts["z_0"]),
        z_1 = FT(opts["z_1"]),
        z_2 = FT(opts["z_2"]),
        rv_0 = FT(opts["rv_0"]),
        rv_1 = FT(opts["rv_1"]),
        rv_2 = FT(opts["rv_2"]),
        tht_0 = FT(opts["tht_0"]),
        tht_1 = FT(opts["tht_1"]),
        tht_2 = FT(opts["tht_2"]),
        precip_sources = Int(opts["precip_sources"]),
        precip_sinks = Int(opts["precip_sinks"]),
        qtot_flux_correction = Int(opts["qtot_flux_correction"]),
        prescribed_Nd = FT(opts["prescribed_Nd"]),
        open_system_activation = Int(opts["open_system_activation"]),
        local_activation = Int(opts["local_activation"]),
        r_dry = FT(opts["r_dry"]),
        std_dry = FT(opts["std_dry"]),
        κ = FT(opts["kappa"]),
    )
    if opts["rain_formation_scheme_choice"] == "SB2006NL"
        toml_dict.data["SB2006_raindrops_terminal_velocity_coeff_aR"]["value"] = FT(6.0)
        toml_dict.data["SB2006_raindrops_terminal_velocity_coeff_bR"]["value"] = FT(9.76)
        toml_dict.data["SB2006_raindrops_terminal_velocity_coeff_cR"]["value"] = FT(1490.0)
    end
    # Create Thermodynamics.jl and KinematicDriver model parameters
    # (some of the CloudMicrophysics.jl parameters structs are created later based on model choices)
    common_params = create_common_parameters(toml_dict)
    kid_params = create_kid_parameters(toml_dict)
    thermo_params = create_thermodynamics_parameters(toml_dict)
    air_params = CMP.AirProperties(toml_dict)
    activation_params = CMP.AerosolActivationParameters(toml_dict)

    moisture = CO.get_moisture_type(moisture_choice, toml_dict)
    precip = CO.get_precipitation_type(
        precipitation_choice, toml_dict;
        rain_formation_choice, sedimentation_choice,
        boundary = opts["p3_boundary_condition"],
    )

    @info "Initialising"
    # Initialize the timestepping struct
    TS = CO.TimeStepping(FT(opts["dt"]), FT(opts["dt_output"]), FT(opts["t_end"]))

    # Create the coordinates
    space, face_space = K1D.make_function_space(FT, FT(opts["z_min"]), FT(opts["z_max"]), opts["n_elem"])
    coord = CC.Fields.coordinate_field(space)
    face_coord = CC.Fields.coordinate_field(face_space)

    # Initialize the netcdf output Stats struct
    output_nc = joinpath(path, "Output.nc")
    Stats = CO.NetCDFIO_Stats(output_nc, 1.0, parent(face_coord), parent(coord))

    # Solve the initial value problem for density profile
    ρ_profile = CO.ρ_ivp(FT, kid_params, thermo_params)
    # Create the initial condition profiles
    if precipitation_choice == "CloudyPrecip"
        pdist_types = determine_cloudy_disttypes(opts["num_moments"])
        cloudy_params, cloudy_pdists = create_cloudy_parameters(FT, pdist_types)
        ic_1d = CO.initial_condition_1d.(FT, common_params, kid_params, thermo_params, (ρ_profile,), coord.z)
        init = CO.cloudy_initial_condition.((cloudy_pdists,), ic_1d)
    elseif precipitation_choice == "PrecipitationP3"
        cloudy_params = nothing
        (; ice_start, _q_flux, _N_flux, _F_rim, _F_liq, _ρ_r_init) = precip.p3_boundary_condition
        init =
            CO.p3_initial_condition.(
                FT, kid_params, thermo_params, coord.z;
                _q_init = _q_flux, _N_init = _N_flux, _F_rim = _F_rim, _F_liq = _F_liq,
                _ρ_r = _ρ_r_init, z_top = FT(opts["z_max"]), ice_start = ice_start,
            )
    else
        cloudy_params = nothing
        init = CO.initial_condition_1d.(FT, common_params, kid_params, thermo_params, (ρ_profile,), coord.z)
    end

    # Create aux vector and apply initial condition
    aux = K1D.initialise_aux(
        FT, init, common_params, kid_params, thermo_params, air_params, activation_params,
        TS, Stats, face_space, moisture, precip, cloudy_params,
    )

    # Create state vector and apply initial condition
    Y = CO.initialise_state(moisture, precip, init)

    # Output the initial condition
    CO.simulation_output(aux, 0.0)

    # Define callbacks for output
    callback_io = ODE.DiscreteCallback(CO.condition_io, CO.affect_io!; save_positions = (false, false))
    callback = ODE.CallbackSet(callback_io)

    # Collect all the tendencies into rhs function for ODE solver
    # based on model choices for the solved equations
    ode_rhs! = K1D.make_rhs_function(moisture, precip)

    # Solve the ODE operator
    problem = ODE.ODEProblem(ode_rhs!, Y, (FT(opts["t_ini"]), FT(opts["t_end"])), aux)
    @info "Solving"
    solver = ODE.solve(
        problem, ODE.SSPRK33();
        TS.dt, saveat = TS.dt_io, callback,
        progress = true, progress_message = (dt, u, p, t) -> t,
    )

    # Some basic plots
    opts["plotting_flag"] == true && with_theme(theme_minimal()) do
        @info "Plotting"
        output = joinpath(path, "figures")

        z_centers = vec(CC.Fields.coordinate_field(space))
        plot_final_aux_profiles(z_centers, aux, precip; output)
        if precip isa CO.PrecipitationP3
            plot_animation_p3(z_centers, solver, aux, moisture, precip, K1D, output)
            plot_timeheight_p3(output_nc, precip; output)
        else
            plot_animation(output_nc; output)
            plot_timeheight(output_nc; output, mixed_phase = false)
        end
    end

    nothing
end
