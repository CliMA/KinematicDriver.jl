"""
    NetCDF output
"""

mutable struct NetCDFIO_Stats
    root_grp::NC.NCDataset{Nothing}
    profiles_grp::NC.NCDataset{NC.NCDataset{Nothing}}
    ts_grp::NC.NCDataset{NC.NCDataset{Nothing}}
    output_interval::Float64
    nc_filename::String
    vars::Dict{String, Any} # Hack to avoid https://github.com/Alexander-Barth/NCDatasets.jl/issues/135
end

function NetCDFIO_Stats(nc_filename, output_interval, z_faces, z_centers, 
    profile_fields; ts_fields=())
    FT = Float64

    # Initialize properties with valid type:
    tmp = tempname()
    root_grp = NC.Dataset(tmp, "c")
    NC.defGroup(root_grp, "profiles")
    NC.defGroup(root_grp, "timeseries")
    profiles_grp = root_grp.group["profiles"]
    ts_grp = root_grp.group["timeseries"]
    close(root_grp)

    # Remove the NC file if it exists, in case it accidentally wasn't closed
    isfile(nc_filename) && rm(nc_filename; force = true)

    NC.Dataset(nc_filename, "c") do root_grp
        # Set profile dimensions
        profile_grp = NC.defGroup(root_grp, "profiles")
        NC.defDim(profile_grp, "zf", length(z_faces))
        NC.defDim(profile_grp, "zc", length(z_centers))
        NC.defDim(profile_grp, "t", Inf)
        NC.defVar(profile_grp, "zf", z_faces, ("zf",))
        NC.defVar(profile_grp, "zc", z_centers, ("zc",))
        NC.defVar(profile_grp, "t", Float64, ("t",))

        for profile_field in Set(profile_fields)
            NC.defVar(profile_grp, profile_field, FT, ("zc", "t"))
        end

        reference_grp = NC.defGroup(root_grp, "reference")
        NC.defDim(reference_grp, "zf", length(z_faces))
        NC.defDim(reference_grp, "zc", length(z_centers))
        NC.defVar(reference_grp, "zf", z_faces, ("zf",))
        NC.defVar(reference_grp, "zc", z_centers, ("zc",))

        ts_grp = NC.defGroup(root_grp, "timeseries")
        NC.defDim(ts_grp, "t", Inf)
        NC.defVar(ts_grp, "t", Float64, ("t",))
        if ~isempty(ts_fields)
            for ts_field in ts_fields
                NC.defVar(ts_grp, ts_field, Float64, ("t",))
            end
        end

    end
    vars = Dict{String, Any}()
    return NetCDFIO_Stats(root_grp, profiles_grp, ts_grp, output_interval, nc_filename, vars)
end

function open_files(self::NetCDFIO_Stats)
    self.root_grp = NC.Dataset(self.nc_filename, "a")
    self.profiles_grp = self.root_grp.group["profiles"]
    self.ts_grp = self.root_grp.group["timeseries"]
    vars = self.vars

    # Hack to avoid https://github.com/Alexander-Barth/NCDatasets.jl/issues/135
    vars["profiles"] = Dict{String, Any}()
    for k in keys(self.profiles_grp)
        vars["profiles"][k] = self.profiles_grp[k]
    end
    vars["timeseries"] = Dict{String, Any}()
    for k in keys(self.ts_grp)
        vars["timeseries"][k] = self.ts_grp[k]
    end
end

function close_files(self::NetCDFIO_Stats)
    close(self.root_grp)
end

function write_field(self::NetCDFIO_Stats, var_name::String, data::T, group) where {T <: AbstractArray{Float64, 1}}
    # Hack to avoid https://github.com/Alexander-Barth/NCDatasets.jl/issues/135
    @inbounds self.vars[group][var_name][:, end] = data
    # Ideally, we remove self.vars and use:
    # var = self.profiles_grp[var_name]
    # Not sure why `end` instead of `end+1`, but `end+1` produces garbage output
    # @inbounds var[end, :] = data :: T
end

function write_ts(self::NetCDFIO_Stats, var_name::String, data::Float64)
    # Hack to avoid https://github.com/Alexander-Barth/NCDatasets.jl/issues/135
    @inbounds self.vars["timeseries"][var_name][end] = data::Float64
    # Ideally, we remove self.vars and use:
    # var = self.ts_grp[var_name]
    # @inbounds var[end+1] = data :: Float64
end

function write_simulation_time(self::NetCDFIO_Stats, t::Float64)
    # # Write to profiles group
    profile_t = self.profiles_grp["t"]
    @inbounds profile_t[end + 1] = t::Float64

    # # Write to timeseries group
    ts_t = self.ts_grp["t"]
    @inbounds ts_t[end + 1] = t::Float64
end

function KiD_output(aux, t::Float64)

    # TODO: remove `vars` hack that avoids
    # https://github.com/Alexander-Barth/NCDatasets.jl/issues/135
    # opening/closing files every step should be okay. #removeVarsHack

    UnPack.@unpack Stats, field_outputs, ts_outputs = aux.io_info
    UnPack.@unpack ρ, p, θ_liq_ice = aux.constants
    UnPack.@unpack T, θ_dry, q_tot, q_liq, q_ice = aux.moisture_variables
    UnPack.@unpack q_rai, q_sno = aux.precip_variables

    open_files(Stats)



    write_simulation_time(Stats, t)

    # TODO: there is probably a more elegant way to do this using dictionaries, but I haven't figured it out
    # field_dict = fieldname_to_var()
    # for field_out in field_outputs
    #     write_field(Stats, field_out, vec(eval(Meta.parse(field_dict[field_out]*"()"))), "profiles")
    # end

    if "density" in field_outputs
        write_field(Stats, "density", vec(ρ), "profiles")
    end
    if "temperature" in field_outputs
        write_field(Stats, "temperature", vec(T), "profiles")
    end
    if "pressure" in field_outputs
        write_field(Stats, "pressure", vec(p), "profiles")
    end
    if "theta_liq_ice" in field_outputs
        write_field(Stats, "theta_liq_ice", vec(ρ), "profiles")
    end
    if "theta_dry" in field_outputs
        write_field(Stats, "density", vec(ρ), "profiles")
    end
    if "q_tot" in field_outputs
        write_field(Stats, "q_tot", vec(q_tot), "profiles")
    end
    if "q_liq" in field_outputs
        write_field(Stats, "q_liq", vec(q_liq), "profiles")
    end
    if "q_ice" in field_outputs
        write_field(Stats, "q_ice", vec(q_ice), "profiles")
    end
    if "q_rai" in field_outputs
        write_field(Stats, "q_rai", vec(q_rai), "profiles")
    end
    if "q_sno" in field_outputs
        write_field(Stats, "q_sno", vec(q_sno), "profiles")
    end

    close_files(Stats)
end

# TODO: make a more elegant dictionary solution for selecting netcdf outputs
# function fieldname_to_var()
#     field_dict = Dict([("density", "ρ")]) #, ("temperature", T), ("pressure", p), 
#         # ("theta_liq_ice",θ_liq_ice), ("theta_dry", θ_dry), ("q_tot", q_tot),
#         # ("q_liq", q_liq), ("q_ice", q_ice), ("q_rai", q_rai), ("q_sno", q_sno)])

#     return field_dict
# end
