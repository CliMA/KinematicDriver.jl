"""
    NetCDF output
"""

struct NetCDFIO_Stats
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

    # Remove the NC file if it exists, in case it accidentally wasn't closed
    isfile(nc_filename) && rm(nc_filename; force = true)

    NC.Dataset(nc_filename, "c") do ds
        # Set profile dimensions
        profiles_grp = NC.defGroup(ds, "profiles")
        NC.defDim(profiles_grp, "zf", length(z_faces))
        NC.defDim(profiles_grp, "zc", length(z_centers))
        NC.defDim(profiles_grp, "t", Inf)
        NC.defVar(profiles_grp, "zf", z_faces, ("zf",))
        NC.defVar(profiles_grp, "zc", z_centers, ("zc",))
        NC.defVar(profiles_grp, "t", FT, ("t",))

        for var in output_profiles
            NC.defVar(profiles_grp, var, FT, ("zc", "t"))
        end

        reference_grp = NC.defGroup(ds, "reference")
        NC.defDim(reference_grp, "zf", length(z_faces))
        NC.defDim(reference_grp, "zc", length(z_centers))
        NC.defVar(reference_grp, "zf", z_faces, ("zf",))
        NC.defVar(reference_grp, "zc", z_centers, ("zc",))

        timeseries_grp = NC.defGroup(ds, "timeseries")
        NC.defDim(timeseries_grp, "t", Inf)
        NC.defVar(timeseries_grp, "t", FT, ("t",))

        #output_timeseries = ["TODO - add timeseries output"]
        #for var in output_timeseries
        #    NC.defVar(timeseries_grp, var, FT, ("t",))
        #end
    end
    return NetCDFIO_Stats(output_interval, nc_filename)
end

function KiD_output(aux, t::Float64)

    (; Stats) = aux

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

    NC.Dataset(Stats.nc_filename, "a") do ds

        # write simulation time to profiles and timeseries
        profile_t = ds.group["profiles"]["t"]
        @inbounds profile_t[end + 1] = t::Float64
        ts_t = ds.group["timeseries"]["t"]
        @inbounds ts_t[end + 1] = t::Float64

        # write profiles
        for field in (aux.constants, aux.moisture_variables, aux.precip_variables)
            for var in keys(name_map)
                if var in propertynames(field)
                    prop = getproperty(field, var)
                    # _parentname = nothing to avoid https://github.com/Alexander-Barth/NCDatasets.jl/issues/135
                    ncdf_var = NC.cfvariable(ds.group["profiles"], name_map[var], _parentname = nothing)
                    @inbounds ncdf_var[:, end + 1] = vec(prop)
                end
            end
        end
        # close the NC dataset
    end
end

