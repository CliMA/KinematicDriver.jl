import NCDatasets as NC

function generate_fake_pysdm_data(dirs; n_files = 50, z_max = 3000.0, n_z = 60, t_max = 240.0, n_t = 40)

    for dir in dirs
        if ~isdir(dir)
            mkdir(dir)
        end
        for i in 1:n_files
            filename = dir * "output_" * string(i) * ".nc"
            NC.Dataset(filename, "c") do ds
                # Define dimensions
                NC.defDim(ds, "height", 60)
                NC.defDim(ds, "time", 40)
                # Define variables:
                ncvar = NC.defVar(ds, "height", Float64, ("height",))
                ncvar[:] = collect(range(0.0, 3000.0, 60))
                ncvar = NC.defVar(ds, "time", Float64, ("time",))
                ncvar[:] = collect(range(0.0, 240.0, 40))
                ncvar = NC.defVar(ds, "qv", Float64, ("time", "height"))
                ncvar[:, :] = rand(40, 60)
                ncvar = NC.defVar(ds, "ql", Float64, ("time", "height"))
                ncvar[:, :] = rand(40, 60)
                ncvar = NC.defVar(ds, "qr", Float64, ("time", "height"))
                ncvar[:, :] = rand(40, 60)
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
