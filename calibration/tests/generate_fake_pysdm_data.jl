import NCDatasets as NC

function generate_fake_pysdm_data(dirs; n_files = 50, z_max = 3000.0, n_z = 60, t_max = 240.0, n_t = 41)

    dz = z_max / n_z
    for dir in dirs
        if ~isdir(dir)
            mkdir(dir)
        end
        for i in 1:n_files
            filename = dir * "output_" * string(i) * ".nc"
            NC.Dataset(filename, "c") do ds
                # Define dimensions
                NC.defDim(ds, "height", n_z)
                NC.defDim(ds, "time", n_t)
                # Define variables:
                ncvar = NC.defVar(ds, "height", Float64, ("height",))
                ncvar[:] = collect(range(dz / 2, z_max - dz / 2, n_z))
                ncvar = NC.defVar(ds, "time", Float64, ("time",))
                ncvar[:] = collect(range(0.0, t_max, n_t))
                ncvar = NC.defVar(ds, "qv", Float64, ("time", "height"))
                ncvar[:, :] = rand(n_t, n_z)
                ncvar = NC.defVar(ds, "qc", Float64, ("time", "height"))
                ncvar[:, :] = rand(n_t, n_z)
                ncvar = NC.defVar(ds, "qr", Float64, ("time", "height"))
                ncvar[:, :] = rand(n_t, n_z)
                close(ds)
            end
        end
    end
end

function rm_fake_pysdm_data(dirs)
    for dir in dirs
        rm(dir, recursive = true, force = true)
    end
end
