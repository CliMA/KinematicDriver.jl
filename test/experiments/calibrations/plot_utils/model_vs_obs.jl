using Plots

import Kinematic1D
const KID = Kinematic1D

import CloudMicrophysics

const CM = CloudMicrophysics
const CMT = CM.CommonTypes
const CM1 = CM.Microphysics1M
const CM2 = CM.Microphysics2M

include("./config.jl")

"""
Define required functions

"""

function non_normalize!(vec, config)
    n_c = length(config["observations"]["cases"])
    n_v = length(config["observations"]["data_names"])
    n_h = config["model"]["n_elem"]
    n_t = length(config["model"]["t_calib"])
    n_vht = n_v * n_h * n_t

    for i in 1:n_c
        v_ = vec[((i - 1) * n_vht + 1):(i * n_vht)]
        m_ = reshape(v_, n_h * n_v, n_t)
        for j in 1:n_v
            m_[((j - 1) * n_h + 1):(j * n_h), :] =
                m_[((j - 1) * n_h + 1):(j * n_h), :] * config["observations"]["ynorm"][i, j]
        end
        v_ = reshape(m_, :, 1)
        vec[((i - 1) * n_vht + 1):(i * n_vht)] = v_
    end
end

function find_max(vec, config)
    n_c = length(config["observations"]["cases"])
    n_h = config["model"]["n_elem"]
    n_t = length(config["model"]["t_calib"])
    n_ht = n_h * n_t
    max_val = zeros(n_c)
    max_t = zeros(n_c)

    for i in 1:n_c
        v_ = vec[((i - 1) * n_ht + 1):(i * n_ht)]
        m_ = reshape(v_, n_h, n_t)
        max_val[i] = maximum(m_)
        for j in 1:n_h
            for k in 1:n_t
                if m_[j, k] == max_val[i]
                    max_t[i] = (k - 1) * 2.5
                    break
                end
            end
        end
    end
    return max_val, max_t
end

function lwp_rwp(vec, config; isref = false)
    n_c = length(config["observations"]["cases"])
    n_v = length(config["observations"]["data_names"])
    n_h = config["model"]["n_elem"]
    n_t = length(config["model"]["t_calib"])
    n_vht = n_v * n_h * n_t
    dz = (config["model"]["z_max"] - config["model"]["z_min"]) / n_h
    lwp = zeros(n_t, n_c)
    rwp = zeros(n_t, n_c)
    rr = zeros(n_t, n_c)

    for i in 1:n_c
        v_ = vec[((i - 1) * n_vht + 1):(i * n_vht)]
        m_ = reshape(v_, n_v * n_h, n_t)
        ρ = m_[1:n_h, :]
        ql = m_[(n_h + 1):(2 * n_h), :]
        qr = m_[(2 * n_h + 1):(3 * n_h), :]
        lwp[:, i] = sum(ql .* ρ, dims = 1) .* dz
        rwp[:, i] = sum(qr .* ρ, dims = 1) .* dz

        ρ0 = ρ[1, :]
        qr0 = qr[1, :]
        params_calib = KID.create_param_dict(u_values, u_names)
        params = KID.create_parameter_set(Float64, config["model"], params_calib)
        microphys_params = KID.Parameters.microphysics_params(params)
        vt0 = CM1.terminal_velocity(microphys_params, CMT.RainType(), ρ, qr0)
        if isref
            vt0 = m_[3 * n_h + 1, :]
        end
        rr[:, i] = qr0 .* ρ0 .* vt0 .* 3600
    end
    return lwp, rwp, rr
end

"""
Read observation and produce model data

"""

data_save_directory = KID.make_output_directories()

config = get_config()
u_names = collect(keys(config["prior"]["parameters"]))
u_values = collect([v.mean for v in values(config["prior"]["parameters"])])

truth = KID.get_obs!(config)
G = KID.run_dyn_model(u_values, u_names, config)

vec_truth = truth.mean * 1000.0
vec_G = G[:] * 1000.0
vec_truth = truth.mean
vec_G = G[:]
non_normalize!(vec_truth, config)
non_normalize!(vec_G, config)

# Rain maximum values and times
max_val_G, max_t_G = find_max(vec_G, config)
max_val_truth, max_t_truth = find_max(vec_truth, config)
max_t_G[25] = 2 * max_t_G[21] - max_t_G[17]
max_t_truth[25] = 2 * max_t_truth[21] - max_t_truth[17]

# LWP, RWP, and rain rate
lwp_G, rwp_G, rr_G = lwp_rwp(vec_G, config)
lwp_truth, rwp_truth, rr_truth = lwp_rwp(vec_truth, config; isref = true)

"""
Plot contours of rain maximum values and times 

"""

Nd = [10, 20, 50, 100, 200, 500, 1000]
rw = [1, 2, 3, 4]

fig = contourf(
    Nd,
    rw,
    reshape(max_val_G, 4, 7),
    levels = 20,
    xscale = :log,
    title = "rai_max value [q/kg]",
    xlabel = "N_d [1/cm³]",
    ylabel = "(ρw)₀ [kg / m² / s]",
    size = (800, 500),
    left_margin = 3Plots.mm,
    bottom_margin = 2Plots.mm,
)
Plots.pdf(fig, "max_val_model.pdf")

fig = contourf(
    Nd,
    rw,
    reshape(max_t_G, 4, 7),
    levels = 20,
    xscale = :log,
    title = "rai_max time [min]",
    xlabel = "N_d [1/cm³]",
    ylabel = "(ρw)₀ [kg / m² / s]",
    size = (800, 500),
    left_margin = 3Plots.mm,
    bottom_margin = 2Plots.mm,
)
Plots.pdf(fig, "max_time_model.pdf")

fig = contourf(
    Nd,
    rw,
    reshape(max_val_truth, 4, 7),
    levels = 20,
    xscale = :log,
    title = "rai_max value [q/kg]",
    xlabel = "N_d [1/cm³]",
    ylabel = "(ρw)₀ [kg / m² / s]",
    size = (800, 500),
    left_margin = 3Plots.mm,
    bottom_margin = 2Plots.mm,
)
Plots.pdf(fig, "max_val_truth.pdf")

fig = contourf(
    Nd,
    rw,
    reshape(max_t_truth, 4, 7),
    levels = 20,
    xscale = :log,
    title = "rai_max time [min]",
    xlabel = "N_d [1/cm³]",
    ylabel = "(ρw)₀ [kg / m² / s]",
    size = (800, 500),
    left_margin = 3Plots.mm,
    bottom_margin = 2Plots.mm,
)
Plots.pdf(fig, "max_time_truth.pdf")


"""
Plot LWP, RWP, and rain rates

"""

time = config["model"]["t_calib"] / 60

plot(
    time,
    lwp_G[:, 4],
    linewidth = 2,
    linecolor = "blue",
    xlabel = "time [min]",
    ylabel = "LWP [kg/m²]",
    label = "(ρw)₀ = 1 kg/m²/s",
)
plot!(time, lwp_G[:, 5], linewidth = 2, linecolor = "red", label = "(ρw)₀ = 2 kg/m²/s")
plot!(time, lwp_G[:, 6], linewidth = 2, linecolor = "green", label = "(ρw)₀ = 3 kg/m²/s")
plot!(time, lwp_G[:, 7], linewidth = 2, linecolor = "magenta", label = "(ρw)₀ = 4 kg/m²/s")
plot!(time, lwp_truth[:, 4], linestyle = :dash, linewidth = 2, linecolor = "blue", label = false)
plot!(time, lwp_truth[:, 5], linestyle = :dash, linewidth = 2, linecolor = "red", label = false)
plot!(time, lwp_truth[:, 6], linestyle = :dash, linewidth = 2, linecolor = "green", label = false)
fig = plot!(
    time,
    lwp_truth[:, 7],
    linestyle = :dash,
    linewidth = 2,
    linecolor = "magenta",
    label = false,
    title = "Nd=100 1/cm³",
)
Plots.pdf(fig, "LWP-rws_Nd=100.pdf")

plot(
    time,
    lwp_G[:, 1],
    linewidth = 2,
    linecolor = "blue",
    xlabel = "time [min]",
    ylabel = "LWP [kg/m²]",
    label = "Nd = 10 1/cm³",
)
plot!(time, lwp_G[:, 2], linewidth = 2, linecolor = "red", label = "Nd = 20 1/cm³")
plot!(time, lwp_G[:, 3], linewidth = 2, linecolor = "yellow", label = "Nd = 50 1/cm³")
plot!(time, lwp_G[:, 6], linewidth = 2, linecolor = "cyan", label = "Nd = 100 1/cm³")
plot!(time, lwp_G[:, 8], linewidth = 2, linecolor = "purple", label = "Nd = 200 1/cm³")
plot!(time, lwp_G[:, 9], linewidth = 2, linecolor = "green", label = "Nd = 500 1/cm³")
plot!(time, lwp_G[:, 10], linewidth = 2, linecolor = "magenta", label = "Nd = 1000 1/cm³")
plot!(time, lwp_truth[:, 1], linestyle = :dash, linewidth = 2, linecolor = "blue", label = false)
plot!(time, lwp_truth[:, 2], linestyle = :dash, linewidth = 2, linecolor = "red", label = false)
plot!(time, lwp_truth[:, 3], linestyle = :dash, linewidth = 2, linecolor = "yellow", label = false)
plot!(time, lwp_truth[:, 6], linestyle = :dash, linewidth = 2, linecolor = "cyan", label = false)
plot!(time, lwp_truth[:, 8], linestyle = :dash, linewidth = 2, linecolor = "purple", label = false)
plot!(time, lwp_truth[:, 9], linestyle = :dash, linewidth = 2, linecolor = "green", label = false)
fig = plot!(
    time,
    lwp_truth[:, 10],
    linestyle = :dash,
    linewidth = 2,
    linecolor = "magenta",
    label = false,
    title = "(ρw)₀ = 3 kg/m²/s",
)
Plots.pdf(fig, "LWP-rw=3_Nds.pdf")

plot(
    time,
    rwp_G[:, 4],
    linewidth = 2,
    linecolor = "blue",
    xlabel = "time [min]",
    ylabel = "RWP [kg/m²]",
    label = "(ρw)₀ = 1 kg/m²/s",
)
plot!(time, rwp_G[:, 5], linewidth = 2, linecolor = "red", label = "(ρw)₀ = 2 kg/m²/s")
plot!(time, rwp_G[:, 6], linewidth = 2, linecolor = "green", label = "(ρw)₀ = 3 kg/m²/s")
plot!(time, rwp_G[:, 7], linewidth = 2, linecolor = "magenta", label = "(ρw)₀ = 4 kg/m²/s")
plot!(time, rwp_truth[:, 4], linestyle = :dash, linewidth = 2, linecolor = "blue", label = false)
plot!(time, rwp_truth[:, 5], linestyle = :dash, linewidth = 2, linecolor = "red", label = false)
plot!(time, rwp_truth[:, 6], linestyle = :dash, linewidth = 2, linecolor = "green", label = false)
fig = plot!(
    time,
    rwp_truth[:, 7],
    linestyle = :dash,
    linewidth = 2,
    linecolor = "magenta",
    label = false,
    title = "Nd=100 1/cm³",
)
Plots.pdf(fig, "RWP-rws_Nd=100.pdf")

plot(
    time,
    rwp_G[:, 1],
    linewidth = 2,
    linecolor = "blue",
    xlabel = "time [min]",
    ylabel = "RWP [kg/m²]",
    label = "Nd = 10 1/cm³",
)
plot!(time, rwp_G[:, 2], linewidth = 2, linecolor = "red", label = "Nd = 20 1/cm³")
plot!(time, rwp_G[:, 3], linewidth = 2, linecolor = "yellow", label = "Nd = 50 1/cm³")
plot!(time, rwp_G[:, 6], linewidth = 2, linecolor = "cyan", label = "Nd = 100 1/cm³")
plot!(time, rwp_G[:, 8], linewidth = 2, linecolor = "purple", label = "Nd = 200 1/cm³")
plot!(time, rwp_G[:, 9], linewidth = 2, linecolor = "green", label = "Nd = 500 1/cm³")
plot!(time, rwp_G[:, 10], linewidth = 2, linecolor = "magenta", label = "Nd = 1000 1/cm³")
plot!(time, rwp_truth[:, 1], linestyle = :dash, linewidth = 2, linecolor = "blue", label = false)
plot!(time, rwp_truth[:, 2], linestyle = :dash, linewidth = 2, linecolor = "red", label = false)
plot!(time, rwp_truth[:, 3], linestyle = :dash, linewidth = 2, linecolor = "yellow", label = false)
plot!(time, rwp_truth[:, 6], linestyle = :dash, linewidth = 2, linecolor = "cyan", label = false)
plot!(time, rwp_truth[:, 8], linestyle = :dash, linewidth = 2, linecolor = "purple", label = false)
plot!(time, rwp_truth[:, 9], linestyle = :dash, linewidth = 2, linecolor = "green", label = false)
fig = plot!(
    time,
    rwp_truth[:, 10],
    linestyle = :dash,
    linewidth = 2,
    linecolor = "magenta",
    label = false,
    title = "(ρw)₀ = 3 kg/m²/s",
)
Plots.pdf(fig, "RWP-rw=3_Nds.pdf")

plot(
    time,
    rr_G[:, 4],
    linewidth = 2,
    linecolor = "blue",
    xlabel = "time [min]",
    ylabel = "Rain rate [mm/hr]",
    label = "(ρw)₀ = 1 kg/m²/s",
)
plot!(time, rr_G[:, 5], linewidth = 2, linecolor = "red", label = "(ρw)₀ = 2 kg/m²/s")
plot!(time, rr_G[:, 6], linewidth = 2, linecolor = "green", label = "(ρw)₀ = 3 kg/m²/s")
plot!(time, rr_G[:, 7], linewidth = 2, linecolor = "magenta", label = "(ρw)₀ = 4 kg/m²/s")
plot!(time, rr_truth[:, 4], linestyle = :dash, linewidth = 2, linecolor = "blue", label = false)
plot!(time, rr_truth[:, 5], linestyle = :dash, linewidth = 2, linecolor = "red", label = false)
plot!(time, rr_truth[:, 6], linestyle = :dash, linewidth = 2, linecolor = "green", label = false)
fig = plot!(
    time,
    rr_truth[:, 7],
    linestyle = :dash,
    linewidth = 2,
    linecolor = "magenta",
    label = false,
    title = "Nd=100 1/cm³",
)
Plots.pdf(fig, "RR-rws_Nd=100.pdf")

plot(
    time,
    rr_G[:, 1],
    linewidth = 2,
    linecolor = "blue",
    xlabel = "time [min]",
    ylabel = "Rain rate [mm/hr]",
    label = "Nd = 10 1/cm³",
)
plot!(time, rr_G[:, 2], linewidth = 2, linecolor = "red", label = "Nd = 20 1/cm³")
plot!(time, rr_G[:, 3], linewidth = 2, linecolor = "yellow", label = "Nd = 50 1/cm³")
plot!(time, rr_G[:, 6], linewidth = 2, linecolor = "cyan", label = "Nd = 100 1/cm³")
plot!(time, rr_G[:, 8], linewidth = 2, linecolor = "purple", label = "Nd = 200 1/cm³")
plot!(time, rr_G[:, 9], linewidth = 2, linecolor = "green", label = "Nd = 500 1/cm³")
plot!(time, rr_G[:, 10], linewidth = 2, linecolor = "magenta", label = "Nd = 1000 1/cm³")
plot!(time, rr_truth[:, 1], linestyle = :dash, linewidth = 2, linecolor = "blue", label = false)
plot!(time, rr_truth[:, 2], linestyle = :dash, linewidth = 2, linecolor = "red", label = false)
plot!(time, rr_truth[:, 3], linestyle = :dash, linewidth = 2, linecolor = "yellow", label = false)
plot!(time, rr_truth[:, 6], linestyle = :dash, linewidth = 2, linecolor = "cyan", label = false)
plot!(time, rr_truth[:, 8], linestyle = :dash, linewidth = 2, linecolor = "purple", label = false)
plot!(time, rr_truth[:, 9], linestyle = :dash, linewidth = 2, linecolor = "green", label = false)
fig = plot!(
    time,
    rr_truth[:, 10],
    linestyle = :dash,
    linewidth = 2,
    linecolor = "magenta",
    label = false,
    title = "(ρw)₀ = 3 kg/m²/s",
)
Plots.pdf(fig, "RR-rw=3_Nds.pdf")
