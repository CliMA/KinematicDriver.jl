"""
    NetCDF output
"""

mutable struct NetCDFIO_Stats
    root_grp::NC.NCDataset{Nothing}
    output_interval::Float64
    nc_filename::String
end

function NetCDFIO_Stats(
    nc_filename,
    output_interval,
    z_faces,
    z_centers;
    output_profiles = [
        "density",
        "temperature",
        "pressure",
        "theta_liq_ice",
        "theta_dry",
        "q_tot",
        "q_liq",
        "q_ice",
        "q_rai",
        "q_sno",
    ],
)
    FT = Float64

    # Initialize properties with valid type:
    tmp = tempname()
    root_grp = NC.Dataset(tmp, "c")
    NC.defGroup(root_grp, "profiles")
    NC.defGroup(root_grp, "timeseries")
    close(root_grp)

    # Remove the NC file if it exists, in case it accidentally wasn't closed
    isfile(nc_filename) && rm(nc_filename; force = true)

    NC.Dataset(nc_filename, "c") do root_grp
        # Set profile dimensions
        profiles_grp = NC.defGroup(root_grp, "profiles")
        NC.defDim(profiles_grp, "zf", length(z_faces))
        NC.defDim(profiles_grp, "zc", length(z_centers))
        NC.defDim(profiles_grp, "t", Inf)
        NC.defVar(profiles_grp, "zf", z_faces, ("zf",))
        NC.defVar(profiles_grp, "zc", z_centers, ("zc",))
        NC.defVar(profiles_grp, "t", FT, ("t",))

        for var in output_profiles
            NC.defVar(profiles_grp, var, FT, ("zc", "t"))
        end

        reference_grp = NC.defGroup(root_grp, "reference")
        NC.defDim(reference_grp, "zf", length(z_faces))
        NC.defDim(reference_grp, "zc", length(z_centers))
        NC.defVar(reference_grp, "zf", z_faces, ("zf",))
        NC.defVar(reference_grp, "zc", z_centers, ("zc",))

        timeseries_grp = NC.defGroup(root_grp, "timeseries")
        NC.defDim(timeseries_grp, "t", Inf)
        NC.defVar(timeseries_grp, "t", FT, ("t",))

        #output_timeseries = ["TODO - add timeseries output"]
        #for var in output_timeseries
        #    NC.defVar(timeseries_grp, var, FT, ("t",))
        #end
    end

    return NetCDFIO_Stats(root_grp, output_interval, nc_filename)
end

function open_files(self::NetCDFIO_Stats)
    self.root_grp = NC.Dataset(self.nc_filename, "a")
end

function close_files(self::NetCDFIO_Stats)
    close(self.root_grp)
end

function write_profiles(self::NetCDFIO_Stats, var_name::String, data::T) where {T <: AbstractArray{Float64, 1}}
    # _parentname = nothing is needed to avoid https://github.com/Alexander-Barth/NCDatasets.jl/issues/135
    ncdf_var = NC.cfvariable(self.root_grp.group["profiles"], var_name, _parentname = nothing)
    # Not sure why not end+1
    @inbounds ncdf_var[:, end] = data
end

function write_timeseries(self::NetCDFIO_Stats, var_name::String, data::Float64)
    # _parentname = nothing is needed to avoid https://github.com/Alexander-Barth/NCDatasets.jl/issues/135
    ncdf_var = NC.cfvariable(self.root_grp.group["timeseries"], var_name, _parentname = nothing)
    # Not sure why not end+1
    @inbounds ncdf_var[end] = data
end

function write_simulation_time(self::NetCDFIO_Stats, t::Float64)
    # # Write to profiles group
    profile_t = self.root_grp.group["profiles"]["t"]
    @inbounds profile_t[end + 1] = t::Float64

    # # Write to timeseries group
    ts_t = self.root_grp.group["timeseries"]["t"]
    @inbounds ts_t[end + 1] = t::Float64
end

function KiD_output(aux, t::Float64)

    (; Stats) = aux
    open_files(Stats)

    write_simulation_time(Stats, t)

    name_map = Dict(
        :ρ => "density",
        :T => "temperature",
        :p => "pressure",
        :θ_liq_ice => "theta_liq_ice",
        :θ_dry => "theta_dry",
        :q_tot => "q_tot",
        :q_liq => "q_liq",
        :q_ice => "q_ice",
        :q_rai => "q_rai",
        :q_sno => "q_sno",
    )
    for field in (aux.constants, aux.moisture_variables, aux.precip_variables)
        for var in keys(name_map)
            if var in propertynames(field)
                prop = getproperty(field, var)
                write_profiles(Stats, name_map[var], vec(prop))
            end
        end
    end

    close_files(Stats)
end
