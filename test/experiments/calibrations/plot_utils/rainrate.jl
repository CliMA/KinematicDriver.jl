import CloudMicrophysics
using Plots, Statistics
import Kinematic1D

const CM = CloudMicrophysics
const CMT = CM.CommonTypes
const CM1 = CM.Microphysics1M
const KID = Kinematic1D

include("./config.jl")

# Configs
config = get_config()
config["prior"]["parameters"] = Dict(
    "τ_cond_evap" => (mean = 10.0, var = 2.0, lbound = 5.0, ubound = 100.0),
    "τ_acnv_rai" => (mean = 1000.0, var = 1e6, lbound = 500.0, ubound = Inf),
    "q_liq_threshold" => (mean = 5.0e-4, var = 5e-5, lbound = 0.0, ubound = Inf),
    "χa_rai" => (mean = 1.0, var = 1.0, lbound = 0.0, ubound = Inf),
    "χv_rai" => (mean = 1.0, var = 0.05, lbound = 0.0, ubound = Inf),
    "a_vent_rai" => (mean = 1.5, var = 0.3, lbound = 0.0, ubound = Inf),
    "b_vent_rai" => (mean = 0.53, var = 0.1, lbound = 0.0, ubound = Inf),
    "Δm_rai" => (mean = 0.0, var = 0.3, lbound = -3.0, ubound = Inf),
    "Δa_rai" => (mean = 0.0, var = 0.2, lbound = -2.0, ubound = Inf),
    "Δv_rai" => (mean = 0.0, var = 0.05, lbound = -0.5, ubound = Inf),
)
params_calib_names = collect(keys(config["prior"]["parameters"]))
config["model"] = get_model_config(params_calib_names)
config["model"]["rain_formation_choice"] = "CliMA_1M"
config["model"]["t_calib"] = 0:50.0:3600.0
config["observations"]["data_names"] = ["qlr", "rho", "rain averaged terminal velocity"]
root_dir = "/Users/sajjadazimi/Postdoc/Results/01-PySDM_1D_rain_shaft/data/03-p1000/"
config["observations"]["cases"] = [(w1 = 3.0, p0 = 100000.0, Nd = 100 * 1e6, dir = root_dir * "rhow=3.0_Nd=100/")]
config["observations"]["ynorm"] =
    ones(length(keys(config["observations"]["cases"])), length(config["observations"]["data_names"]))

# Obs
dz = (config["model"]["z_max"] - config["model"]["z_min"]) / config["model"]["n_elem"]
heights =
    collect(range(config["model"]["z_min"] + dz / 2, config["model"]["z_max"] - dz / 2, config["model"]["n_elem"]))
times = collect(config["model"]["t_calib"])
variables = config["observations"]["data_names"]
cases = config["observations"]["cases"]
nt = length(times)
nv = length(variables)
nh = config["model"]["n_elem"]
truth = mean(KID.get_obs_matrix(cases, variables, heights, times), dims = 2)
qlr_0_truth = reshape(truth, nh, :)[1, 1:nv:(nv * nt)]
ρ_0_truth = reshape(truth, nh, :)[1, 2:nv:(nv * nt)]
vt_0_truth = reshape(truth, nh, :)[1, 3:nv:(nv * nt)]
R_truth = ρ_0_truth .* qlr_0_truth .* vt_0_truth .* 3600.0

# Compute model 1: CliMA_1M default
u_names = collect(keys(config["prior"]["parameters"]))
u_values = collect([v.mean for v in values(config["prior"]["parameters"])])
G = KID.run_dyn_model(u_values, u_names, config)

rlr_0_G = reshape(G, nh, :)[1, 1:nv:(nv * nt)] * config["observations"]["ynorm"][1]
ρ_0_G = reshape(G, nh, :)[1, 2:nv:(nv * nt)] * config["observations"]["ynorm"][2]
params_calib = KID.create_param_dict(u_values, u_names)
params = KID.create_parameter_set(Float64, config["model"], params_calib)
microphys_params = KID.Parameters.microphysics_params(params)
Base.broadcastable(ps::CMT.VarNrType) = Ref(ps)
vt_0_G = CM1.terminal_velocity.(microphys_params, CMT.RainType(), ρ_0_G, rlr_0_G)
R_1M_d = ρ_0_G .* rlr_0_G .* vt_0_G .* 3600.0

# Compute model 2: CliMA_1M calibrated
config["prior"]["parameters"] = Dict(
    "τ_cond_evap" => (mean = 19.81, var = 2.5, lbound = 5.0, ubound = 100.0),
    "τ_acnv_rai" => (mean = 1.015e7, var = 10.0, lbound = 100.0, ubound = Inf),
    "q_liq_threshold" => (mean = 3.979e-5, var = 5e-5, lbound = 0.0, ubound = Inf),
    "χa_rai" => (mean = 14.25, var = 0.4, lbound = 0.0, ubound = Inf),
    "χv_rai" => (mean = 0.1137, var = 0.015, lbound = 0.0, ubound = Inf),
    "a_vent_rai" => (mean = 0.001395, var = 0.15, lbound = 0.0, ubound = Inf),
    "b_vent_rai" => (mean = 5.403, var = 0.053, lbound = 0.0, ubound = Inf),
    "Δm_rai" => (mean = 1.273, var = 0.3, lbound = -3.0, ubound = Inf),
    "Δa_rai" => (mean = -0.6628, var = 0.2, lbound = -2.0, ubound = Inf),
    "Δv_rai" => (mean = 2.089, var = 0.05, lbound = -0.5, ubound = Inf),
)
u_names = collect(keys(config["prior"]["parameters"]))
u_values = collect([v.mean for v in values(config["prior"]["parameters"])])
G = KID.run_dyn_model(u_values, u_names, config)

rlr_0_G = reshape(G, nh, :)[1, 1:nv:(nv * nt)] * config["observations"]["ynorm"][1]
ρ_0_G = reshape(G, nh, :)[1, 2:nv:(nv * nt)] * config["observations"]["ynorm"][2]
params_calib = KID.create_param_dict(u_values, u_names)
params = KID.create_parameter_set(Float64, config["model"], params_calib)
microphys_params = KID.Parameters.microphysics_params(params)
Base.broadcastable(ps::CMT.VarNrType) = Ref(ps)
vt_0_G = CM1.terminal_velocity.(microphys_params, CMT.RainType(), ρ_0_G, rlr_0_G)
R_1M_c = ρ_0_G .* rlr_0_G .* vt_0_G .* 3600.0

# Plot
t = 0:50:3600
plot(t ./ 60, R_truth, linewidth = 3, label = "PySDM")
plot!(t ./ 60, R_1M_c, linewidth = 3, label = "1M-calibrated", linestyle = :dashdot)
plot!(t ./ 60, R_1M_d, linewidth = 3, label = "1M-default", linestyle = :dot)
plot!(xlabel = "time [min]")
plot!(ylabel = "rainrate [mm/hr]")
plot!(labelfontsize = 12)
plot!(titlefontsize = 12)
plot!(tickfontsize = 10)
plot!(legendfontsize = 11)
# plot!(legend=:topleft)
fig = plot!()
Plots.pdf(fig, "rainrate_model_vs_pysdm.pdf")
