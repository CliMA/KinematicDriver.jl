"""
    NetCDF output
"""

struct NetCDFIO_Stats
    output_interval::Float64
    nc_filename::String
    output_profiles::Dict{Symbol, String}
    output_precip_sources::Dict{Symbol, String}
    output_activation_sources::Dict{Symbol, String}
end

function NetCDFIO_Stats(
    nc_filename,
    output_interval,
    z_faces,
    z_centers;
    x = nothing,
    output_profiles = Dict(
        :ρ => "density",
        :T => "temperature",
        :p => "pressure",
        :θ_liq_ice => "theta_liq_ice",
        :θ_dry => "theta_dry",
        :q_tot => "q_tot",
        :q_liq => "q_liq",
        :q_ice => "q_ice",
        :q_rim => "q_rim",
        :q_liqonice => "q_liqonice",
        :q_rai => "q_rai",
        :q_sno => "q_sno",
        :q_vap => "q_vap",
        :N_aer => "N_aer",
        :N_liq => "N_liq",
        :N_rai => "N_rai",
        :N_ice => "N_ice",
        :B_rim => "B_rim",
        :ρq_ice => "ρq_ice",
        :ρq_rim => "ρq_rim",
        :ρq_liqonice => "ρq_liqonice",
        :ρq_vap => "ρq_vap",
        :ρq_rai => "ρq_rai",
    ),
    output_precip_sources = Dict(
        :q_liq => "Sq_liq_prc",
        :q_rai => "Sq_rai_prc",
        :N_aer => "SN_aer_prc",
        :N_liq => "SN_liq_prc",
        :N_rai => "SN_rai_prc",
    ),
    output_activation_sources = Dict(
        :N_aer => "SN_aer_act",
        :N_liq => "SN_liq_act",
    ),
)
    FT = Float64

    dim = x === nothing ? 1 : 2

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

        if dim == 2
            NC.defDim(profiles_grp, "x", length(x))
            NC.defVar(profiles_grp, "x", x, ("x",))
        end

        for field in (output_profiles, output_precip_sources, output_activation_sources)
            for var in keys(field)
                if dim == 1
                    NC.defVar(profiles_grp, field[var], FT, ("zc", "t"))
                elseif dim == 2
                    NC.defVar(profiles_grp, field[var], FT, ("x", "zc", "t"))
                else
                    error("Invalid system dimension!!!")
                end
            end
        end

        reference_grp = NC.defGroup(ds, "reference")
        NC.defDim(reference_grp, "zf", length(z_faces))
        NC.defDim(reference_grp, "zc", length(z_centers))
        NC.defVar(reference_grp, "zf", z_faces, ("zf",))
        NC.defVar(reference_grp, "zc", z_centers, ("zc",))

        if dim == 2
            NC.defDim(reference_grp, "x", length(x))
            NC.defVar(reference_grp, "x", x, ("x",))
        end

        timeseries_grp = NC.defGroup(ds, "timeseries")
        NC.defDim(timeseries_grp, "t", Inf)
        NC.defVar(timeseries_grp, "t", FT, ("t",))
        NC.defVar(timeseries_grp, "ql_max", FT, ("t",))
    end
    return NetCDFIO_Stats(
        output_interval,
        nc_filename,
        output_profiles,
        output_precip_sources,
        output_activation_sources,
    )
end

"""
    write_fields(field, mapping, group, dim)

Writes selected variables from a given `field` to a NetCDF group using `NCDatasets.jl`. Only variables
listed in `mapping` and present in `field` are written.

# Arguments
- `field`: A `NamedTuple` tupple containing simulation variables (e.g., `aux.thermo_variables`).
- `mapping`: A `Dict` mapping variable names (symbols) to NetCDF variable names (strings), e.g., `Stats.output_profiles`.
- `group`: The NetCDF group to which variables are written, e.g., `ds.group["profiles"]`.
- `dim`: The dimensionality of the simulation (1 or 2), used to select the correct data slicing.
"""
function write_fields(field, mapping, group, dim)
    for var in keys(mapping)
        if var in propertynames(field)
            prop = getproperty(field, var)
            # Avoid issue https://github.com/Alexander-Barth/NCDatasets.jl/issues/135
            ncdf_var = NC.cfvariable(group, mapping[var], _parentname = nothing)
            if dim == 1
                @inbounds ncdf_var[:, end] = vec(prop)
            elseif dim == 2
                @inbounds ncdf_var[:, :, end] = parent(prop)[:, 1, 1, :]
            else
                error("Invalid system dimension!")
            end
        end
    end
end

function simulation_output(aux, t::FT) where {FT <: Real}

    (; Stats) = aux
    dim = :domain_width in keys(aux) ? 2 : 1

    NC.Dataset(Stats.nc_filename, "a") do ds

        # write simulation time to profiles and timeseries
        profile_t = ds.group["profiles"]["t"]
        @inbounds profile_t[end + 1] = t::FT
        ts_t = ds.group["timeseries"]["t"]
        @inbounds ts_t[end + 1] = t::FT

        # write profiles
        grp = ds.group["profiles"]
        write_fields(aux.thermo_variables, Stats.output_profiles, grp, dim)
        write_fields(aux.microph_variables, Stats.output_profiles, grp, dim)
        write_fields(aux.precip_sources, Stats.output_precip_sources, grp, dim)
        write_fields(aux.activation_sources, Stats.output_activation_sources, grp, dim)

        # write timeseries
        ncdf_var = NC.cfvariable(ds.group["timeseries"], "ql_max", _parentname = nothing)
        @inbounds ncdf_var[end] = maximum(aux.microph_variables.q_liq)

        # close the NC dataset
    end
end

"""
   Interface to ODE callbacks for handling output.
"""

function condition_io(u, t, integrator)
    (; TS, Stats) = integrator.p
    TS.dt_io += TS.dt
    io_flag = false
    if TS.dt_io > Stats.output_interval
        TS.dt_io = 0
        io_flag = true
    end
    return io_flag || t ≈ 0 || t ≈ TS.t_max
end

"""
   Interface to ODE callbacks for handling output.
"""
function affect_io!(integrator)
    simulation_output(integrator.p, integrator.t)

    ODE.u_modified!(integrator, false) # We're legitamately not mutating `u` (the state vector)
end
