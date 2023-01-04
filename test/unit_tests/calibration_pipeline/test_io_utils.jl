@testset "make output directory - save - load" begin
    #setup
    tmp_dir = "/tmp/"
    u_names = ["a1", "a2"]
    u_values = [1.0, 2.0]
    u_ensemble = [1.0 1.1 1.2; 2.0 2.1 2.2]
    res = [u_ensemble, u_ensemble]
    config = get_config()
    tmp_file = pwd() * tmp_dir * "tmp_params.jld2"

    #action
    KID.make_output_directories(tmp_dir)
    KID.save_data(res, u_values, u_names, config, file_name = tmp_file)
    tmp_data = KID.load_data(tmp_file)

    #test
    @test isdir(pwd() * tmp_dir)
    @test tmp_data.result == res
    @test tmp_data.u_names == u_names
    @test tmp_data.u_bests == u_values
    @test tmp_data.config == config

    rm(pwd() * tmp_dir, recursive = true, force = true)
end

@testset "plot correlation map" begin
    #setup
    u_names = ["a", "b"]
    u_cov = [0.5 0.6; 0.1 0.2]
    tmp_filename = "tmp.png"

    #action
    KID.plot_correlation_map(u_names, u_cov; output_filename = tmp_filename)

    #test
    @test isfile(tmp_filename)
    rm(tmp_filename)
end
