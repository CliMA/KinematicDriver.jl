import Combinatorics
import NCDatasets as NC
include("../../Calibration.jl")

config_dir = "../K1D"
data_save_directory = make_output_directories(joinpath(@__DIR__, config_dir, "run_model_output"))
output_file_name = data_save_directory * "grid_search.nc"


parameters = Dict(
    "rain_cross_section_size_relation_coefficient_chia" => (ref = 1.0, range = range(0.5, 5.0, 10)),
    "rain_terminal_velocity_size_relation_coefficient_chiv" => (ref = 1.0, range = range(0.2, 2.0, 10)),
)

include(joinpath(@__DIR__, config_dir, "config.jl"))
config = get_config()
param_names = collect(keys(parameters))
model_settings = get_model_config()
config["model"] = model_settings
param_names = collect(keys(parameters))
truth = get_obs!(config)
ref_stats_list = make_ref_stats_list(truth, config["statistics"], get_numbers_from_config(config)...)
ref_stats = combine_ref_stats(ref_stats_list)

NC.Dataset(output_file_name, "c") do ds
    for (combination_num, (param_name_1, param_name_2)) in enumerate(Combinatorics.combinations(param_names, 2))
        group_name = "$param_name_1.$param_name_2"
        println("Combination ", combination_num, " : ", group_name)

        # compute loss, store in matrix (p1, p2)
        value1 = parameters[param_name_1].range
        value2 = parameters[param_name_2].range
        n_value1 = length(value1)
        n_value2 = length(value2)
        loss_2D_sec = zeros((n_value1, n_value2))
        tmp_value_dict = Dict{String, Float64}([Pair(name, value.ref) for (name, value) in parameters])
        for (i, v1) in enumerate(value1)
            tmp_value_dict[param_name_1] = v1
            for (j, v2) in enumerate(value2)
                tmp_value_dict[param_name_2] = v2
                ϕ_names = collect(keys(tmp_value_dict))
                ϕ_values = collect(values(tmp_value_dict))
                loss_2D_sec[i, j, 1, 1] = compute_loss(ϕ_values, ϕ_names, config, ref_stats)
                println("\t i ", i, ", j ", j, " -->  loss = ", loss_2D_sec[i, j])
            end
        end

        # save output
        NC.defGroup(ds, group_name, attrib = [])
        group_root = ds.group[group_name]
        # Define dimensions: param_name_1, param_name_2, case, ensemble_member
        NC.defDim(group_root, "param_name_1", n_value1)
        NC.defDim(group_root, "param_name_2", n_value2)
        # Define variables: 
        # param_name_1 and param_name_2 values
        ncvar = NC.defVar(group_root, param_name_1, value1, ("param_name_1",))
        ncvar[:] = value1
        ncvar = NC.defVar(group_root, param_name_2, value2, ("param_name_2",))
        ncvar[:] = value2
        ncvar = NC.defVar(group_root, "loss_data", loss_2D_sec, ("param_name_1", "param_name_2"))
        ncvar[:, :] = loss_2D_sec
    end  # end for param combinations loop
end  # NC.Dataset do-block
