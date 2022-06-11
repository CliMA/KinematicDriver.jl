"""
Plotting utilities
"""
ENV["GKSwstype"] = "nul"
using ClimaCorePlots, Plots
Plots.GRBackend()

function plot_anim(solver, varselect, aux, z_centers, path)
    anim = Plots.@animate for u in solver.u
        if varselect == "ρq_tot"
            ρq = parent(u.ρq_tot)
        elseif varselect == "ρq_liq"
            ρq = parent(u.ρq_liq)
        elseif varselect == "ρq_ice"
            ρq = parent(u.ρq_ice)
        elseif varselect == "ρq_rai"
            ρq = parent(u.ρq_rai)
        elseif varselect == "ρq_sno"
            ρq = parent(u.ρq_sno)
        else
            error("variable not defined: $varselect")
        end

        ρ = parent(aux.constants.ρ)
        Plots.plot(ρq ./ ρ, z_centers)
    end
    Plots.mp4(anim, joinpath(path, varselect, ".mp4"), fps=10)
end

function plot_spaghetti(data, z_centers, path, filename)
    Plots.png(Plots.plot(data, z_centers), joinpath(path, filename))
end

function plot_timeheight(solver, varselect, aux, z_plot, path)
    t_plot = (t for t in parent(solver.t))
    var_plot = ()
    for u in solver.u
        if varselect == "q_tot"
            ρq = parent(u.ρq_tot)
        elseif varselect == "q_liq"
            ρq = parent(u.ρq_liq)
        elseif varselect == "q_ice"
            ρq = parent(u.ρq_ice)
        elseif varselect == "q_rai"
            ρq = parent(u.ρq_rai)
        elseif varselect == "q_sno"
            ρq = parent(u.ρq_sno)
        else
            error("variable not defined: $varselect")
        end

        ρ = parent(aux.constants.ρ)
        var_plot.append(parent(ρq ./ ρ))
    end

    Plots.png(Plots.heatmap(t_plot, z_plot, var_plot), joinpath(path, varselect, ".png"))
end
# TODO: timeheight for auxiliary variables (?)