"""
    NetCDF output
"""

struct NetCDFIO_Stats
    output_interval::Float64
    nc_filename::String
    output_profiles::Dict{Symbol, String}
end

function NetCDFIO_Stats(
    nc_filename,
    output_interval,
    x,
    z_faces,
    z_centers;
    output_profiles = Dict(
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
    ),
)
    FT = Float64

    # Remove the NC file if it exists, in case it accidentally wasn't closed
    isfile(nc_filename) && rm(nc_filename; force = true)

    NC.Dataset(nc_filename, "c") do ds
        # Set profile dimensions
        profiles_grp = NC.defGroup(ds, "profiles")
        NC.defDim(profiles_grp, "x", length(x))
        NC.defDim(profiles_grp, "zf", length(z_faces))
        NC.defDim(profiles_grp, "zc", length(z_centers))
        NC.defDim(profiles_grp, "t", Inf)
        NC.defVar(profiles_grp, "x", x, ("x",))
        NC.defVar(profiles_grp, "zf", z_faces, ("zf",))
        NC.defVar(profiles_grp, "zc", z_centers, ("zc",))
        NC.defVar(profiles_grp, "t", FT, ("t",))

        for var in keys(output_profiles)
            NC.defVar(profiles_grp, output_profiles[var], FT, ("x", "zc", "t"))
        end

        reference_grp = NC.defGroup(ds, "reference")
        NC.defDim(reference_grp, "x", length(x))
        NC.defDim(reference_grp, "zf", length(z_faces))
        NC.defDim(reference_grp, "zc", length(z_centers))
        NC.defVar(reference_grp, "x", x, ("x",))
        NC.defVar(reference_grp, "zf", z_faces, ("zf",))
        NC.defVar(reference_grp, "zc", z_centers, ("zc",))

        timeseries_grp = NC.defGroup(ds, "timeseries")
        NC.defDim(timeseries_grp, "t", Inf)
        NC.defVar(timeseries_grp, "t", FT, ("t",))
        NC.defVar(timeseries_grp, "ql_max", FT, ("t",))
    end
    return NetCDFIO_Stats(output_interval, nc_filename, output_profiles)
end

function Kinematic2D_output(aux, t::FT) where {FT <: Real}

    (; Stats) = aux

    NC.Dataset(Stats.nc_filename, "a") do ds

        # write simulation time to profiles and timeseries
        profile_t = ds.group["profiles"]["t"]
        @inbounds profile_t[end + 1] = t::FT
        ts_t = ds.group["timeseries"]["t"]
        @inbounds ts_t[end + 1] = t::FT

        # write profiles
        for field in (aux.moisture_variables, aux.precip_variables)
            for var in keys(Stats.output_profiles)
                if var in propertynames(field)
                    prop = getproperty(field, var)
                    # _parentname = nothing to avoid https://github.com/Alexander-Barth/NCDatasets.jl/issues/135
                    ncdf_var = NC.cfvariable(ds.group["profiles"], Stats.output_profiles[var], _parentname = nothing)
                    @inbounds ncdf_var[:, :, end] = parent(prop)[:, 1, 1, :]
                end
            end
        end
        # write timeseries
        ncdf_var = NC.cfvariable(ds.group["timeseries"], "ql_max", _parentname = nothing)
        @inbounds ncdf_var[end] = maximum(aux.moisture_variables.q_liq)

        # close the NC dataset
    end
end

"""
   Interface to ODE callbacks for handling output.
"""
function affect_io!(integrator)
    Kinematic2D_output(integrator.p, integrator.t)

    ODE.u_modified!(integrator, false) # We're legitamately not mutating `u` (the state vector)
end
