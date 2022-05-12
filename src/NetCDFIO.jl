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

function NetCDFIO_Stats(nc_filename, output_interval, z_faces, z_centers)
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

        # TODO - define output variables based on the model that is being run?
        NC.defVar(profile_grp, "density", FT, ("zc","t"))
        NC.defVar(profile_grp, "temperature", FT, ("zc","t"))
        NC.defVar(profile_grp, "pressure", FT, ("zc","t"))

        NC.defVar(profile_grp, "theta_liq_ice", FT, ("zc","t"))
        NC.defVar(profile_grp, "theta_dry", FT, ("zc","t"))

        NC.defVar(profile_grp, "q_tot", FT, ("zc","t"))
        NC.defVar(profile_grp, "q_liq", FT, ("zc","t"))
        NC.defVar(profile_grp, "q_ice", FT, ("zc","t"))
        NC.defVar(profile_grp, "q_rai", FT, ("zc","t"))
        NC.defVar(profile_grp, "q_sno", FT, ("zc","t"))

        reference_grp = NC.defGroup(root_grp, "reference")
        NC.defDim(reference_grp, "zf", length(z_faces))
        NC.defDim(reference_grp, "zc", length(z_centers))
        NC.defVar(reference_grp, "zf", z_faces, ("zf",))
        NC.defVar(reference_grp, "zc", z_centers, ("zc",))

        ts_grp = NC.defGroup(root_grp, "timeseries")
        NC.defDim(ts_grp, "t", Inf)
        NC.defVar(ts_grp, "t", Float64, ("t",))
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

    UnPack.@unpack Stats, ρ, T, θ_dry, p, q_tot, q_liq, q_ice, q_rai, q_sno, θ_liq_ice = aux

    open_files(Stats)

    write_simulation_time(Stats, t)
    write_field(Stats, "density", vec(ρ), "profiles")
    write_field(Stats, "temperature", vec(T), "profiles")
    write_field(Stats, "pressure", vec(p), "profiles")

    write_field(Stats, "theta_liq_ice", vec(θ_liq_ice), "profiles")
    write_field(Stats, "theta_dry", vec(θ_dry), "profiles")

    write_field(Stats, "q_tot", vec(q_tot), "profiles")
    write_field(Stats, "q_liq", vec(q_liq), "profiles")
    write_field(Stats, "q_ice", vec(q_ice), "profiles")
    write_field(Stats, "q_rai", vec(q_rai), "profiles")
    write_field(Stats, "q_sno", vec(q_sno), "profiles")

    close_files(Stats)
end
