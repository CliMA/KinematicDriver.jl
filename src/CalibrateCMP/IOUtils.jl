function make_output_directories(; dir::String = "/output/")
    homedir = pwd()
    data_save_directory = homedir * dir
    if ~isdir(data_save_directory)
        mkdir(data_save_directory)
    end
    return data_save_directory
end

function save_data(res, u_bests, u_names, config; file_name = "parameters.jld2")
    @save file_name res u_bests u_names config
end

function load_data(file_name::String = "parameters.jld2")
    @load file_name res u_bests u_names config
    return (result = res, u_bests = u_bests, u_names = u_names, config = config)
end

function print_results(u, u_names)
    println("Optimization results:")
    for (i, name) in enumerate(u_names)
        println(name, " = ", round(u[i], sigdigits = 4))
    end
end

function ensemble_convergence_single_panel(
    u::Vector{Matrix{FT}},
    sol;
    file_name::String = "ensemble_convergence.gif",
    sol_label = "true_parameters",
    params_index = [1, 2],
    axis_labels = ["u₁", "u₂"],
) where {FT <: Real}
    n_iter_p1 = length(u)
    limits = get_limits(u)

    anim_conv = @animate for i in 1:n_iter_p1
        u_i = u[i]

        plot(
            [sol[params_index[1]]],
            [sol[params_index[2]]],
            seriestype = :scatter,
            markershape = :star5,
            markersize = 11,
            markercolor = :red,
            label = sol_label,
        )

        plot!(
            u_i[params_index[1], :],
            u_i[params_index[2], :],
            seriestype = :scatter,
            xlims = limits[params_index[1]],
            ylims = limits[params_index[2]],
            xlabel = axis_labels[1],
            ylabel = axis_labels[2],
            markersize = 5,
            markeralpha = 0.6,
            markercolor = :blue,
            label = "particles",
            title = "EKI iteration = " * string(i),
        )
    end

    gif(anim_conv, file_name, fps = 2)
end

function ensemble_convergence(
    u::Vector{Matrix{FT}},
    priors,
    config::Dict;
    file_name::String = "ensemble_convergence.gif",
) where {FT <: Real}
    θ = params_validation(config, priors).θ
    ensemble_convergence(u, θ, file_name = file_name)
end

function ensemble_convergence(
    u::Vector{Matrix{FT}},
    sol::Array{Float64, 1};
    file_name::String = "ensemble_convergence.gif",
) where {FT <: Real}

    n_iter_p1 = length(u)
    n_par = size(u[1])[1]
    n_panel = Int(ceil(n_par / 2))
    limits = get_limits(u)

    n_panel_r = Int(floor(sqrt(n_panel)))
    n_panel_c = Int(ceil(sqrt(n_panel)))
    if (n_panel_r * n_panel_c < n_panel)
        n_panel_r = n_panel_r + 1
    end
    layout = (n_panel_r, n_panel_c)
    font_size = n_panel_c < 4 ? Int(20 - 4 * n_panel_c) : 6

    anim_conv = @animate for i in 1:n_iter_p1
        u_i = u[i]
        p = Array{Plots.Plot}(undef, n_panel)

        for j in 1:n_panel
            index_1 = 2 * j - 1
            index_2 = (j < n_panel || iseven(n_par)) ? 2 * j : 1

            p[j] = plot(
                [sol[index_1]],
                [sol[index_2]],
                seriestype = :scatter,
                markershape = :star5,
                markersize = 11,
                markercolor = :red,
            )

            plot!(
                p[j],
                u_i[index_1, :],
                u_i[index_2, :],
                seriestype = :scatter,
                xlims = limits[index_1],
                ylims = limits[index_2],
                xlabel = "u" * string(index_1),
                ylabel = "u" * string(index_2),
                markersize = 5,
                markeralpha = 0.6,
                markercolor = :blue,
                ticks = [],
            )
        end

        plot(
            p...,
            layout = layout,
            labelfontsize = font_size,
            title = "iteration = " * string(i),
            titlefontsize = font_size,
            legend = false,
            framestyle = :box,
            widen = false,
            left_margin = 5Plots.mm,
        )

    end

    gif(anim_conv, file_name, fps = 2)
end

function read_pysdm_data(file)

    @assert endswith(file, ".nc")

    ds = NC.Dataset(file, "r")

    # iterate over all attributes
    attrs = Dict()
    for (attrname, attrval) in ds.attrib
        attrs[attrname] = attrval
    end

    # iterate over all variables
    vars = Dict()
    for (varname, var) in ds
        vars[varname] = reshape(var[:], size(var))
    end

    close(ds)

    data = (attributes = attrs, variables = vars)
    return data
end

function compare_model_and_obs_contours(
    G::Vector{FT},
    obs::Vector{FT},
    config::Dict;
    variables = config["observations"]["data_names"],
    levels = 10,
    linewidth = 0.5,
    fontsize = 10,
    G_title = "model",
    obs_title = "observation",
    path = "output",
    file_base = "model_vs_obs_contours",
) where {FT <: Float64}

    (; n_cases, n_heights, n_times) = get_numbers_from_config(config)

    _dt = (config["model"]["t_end"] - config["model"]["t_ini"]) / length(config["model"]["t_calib"])
    _times_f::Array{FT} = collect(config["model"]["t_calib"])
    _times_c::Array{FT} = collect(
        range(
            config["model"]["t_ini"] + _dt / 2,
            config["model"]["t_end"] - _dt / 2,
            length(config["model"]["t_calib"]) - 1,
        ),
    )
    _times::Array{FT} = config["model"]["filter"]["apply"] == true ? _times_c : _times_f
    _variables_full = config["observations"]["data_names"]
    @assert issubset(Set(variables), Set(_variables_full))
    n_variables = length(variables)
    n_variables_full = length(_variables_full)
    @assert n_variables_full == length(n_heights)
    n_single_case = sum(n_heights) * n_times

    layout = (2, n_variables)
    titles = vcat([v * "_" * G_title for v in variables], [v * "_" * obs_title for v in variables])
    
    for case_num in 1:n_cases
        _v_model = get_case_i_vec(G, case_num, n_single_case)
        _v_obs = get_case_i_vec(obs, case_num, n_single_case)
        _fields_model = get_single_case_fields(_v_model, n_heights, n_times)
        _fields_obs = get_single_case_fields(_v_obs, n_heights, n_times)

        p = Array{Plots.Plot}(undef, n_variables * 2)

        for (i, _data) in enumerate([_fields_model, _fields_obs])
            k = 1
            for j in 1:n_variables_full
                if !(_variables_full[j] in variables)
                    continue
                end

                # find z levels for variable j
                _dz = (config["model"]["z_max"] - config["model"]["z_min"]) / n_heights[j]
                _heights::Array{FT} =
                    collect(range(config["model"]["z_min"] + _dz / 2, config["model"]["z_max"] - _dz / 2, n_heights[j]))
                _heights_f::Array{FT} = collect(range(config["model"]["z_min"], config["model"]["z_max"], n_heights[j] + 1))
                
                # find data for variable j
                _data_var = _data[j]
                if n_heights[j] == 1
                    _f_1d = LinearInterpolation(_times, _data_var[:], extrapolation_bc = Line())
                    _f(z, t) = _f_1d(t)
                else
                    _f = LinearInterpolation((_heights, _times), _data_var, extrapolation_bc = Line())
                end
                _data_f = [_f(z, t) for z in _heights_f, t in _times_f]
                _data_c = [_f(z, t) for z in _heights, t in _times_c]

                if config["model"]["filter"]["apply"] == true
                    p[(i - 1) * n_variables + k] = heatmap(
                        _times_c ./ 60,
                        _heights ./ 1000,
                        _data_c,
                        xlabel = "time [min]",
                        ylabel = "height [km]",
                        title = titles[(i - 1) * n_variables + k],
                    )
                else
                    p[(i - 1) * n_variables + k] = contourf(
                        _times_f ./ 60,
                        _heights_f ./ 1000,
                        _data_f,
                        levels = levels,
                        linewidth = linewidth,
                        xlabel = "time [min]",
                        ylabel = "height [km]",
                        title = titles[(i - 1) * n_variables + k],
                    )
                end
                k = k + 1
            end
        end
        
        fig = plot(
            p...,
            layout = layout,
            labelfontsize = fontsize,
            titlefontsize = fontsize,
            size = (n_variables * 675, 750),
            left_margin = n_variables * 3Plots.mm,
            right_margin = 0Plots.mm,
            bottom_margin = n_variables * 2Plots.mm,
        )
        path = joinpath(@__DIR__, path)
        mkpath(path)
        file_name = file_base * "_case_" * string(case_num) * ".png"
        Plots.png(fig, joinpath(path, file_name))
    end
end

"""
    plot_correlation_map(file_uki; output_filename)

plot correlation map and save as a png file. 

# Inputs:
- `u_names` :: parameter names
- `u_cov` :: parameter covariance matrix
- `output_filename` :: name of the output png file
"""
function plot_correlation_map(
    u_names::Vector{String},
    u_cov::Matrix{FT};
    output_filename::String = "correlation_map.png",
) where {FT <: Real}
    num = size(u_cov)[1]
    cor_map = zeros(num, num)
    for i in 1:num
        for j in num:-1:1
            if j > num - i
                cor_map[i, j] = NaN
            else
                cor_map[i, j] = u_cov[i, num - j + 1] / sqrt(u_cov[i, i]) / sqrt(u_cov[num - j + 1, num - j + 1])
            end
        end
    end

    fig = heatmap(
        cor_map,
        xticks = (collect(1:num), collect(u_names)[end:-1:2]),
        yticks = (collect(1:num), collect(u_names)[1:(end - 1)]),
        size = (800, 600),
    )
    Plots.png(fig, output_filename)
end
