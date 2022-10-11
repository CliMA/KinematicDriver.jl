@testset "make output directory - save - load" begin
    #setup
    tmp_dir = "/tmp/"
    u_names = ["a1", "a2"]
    u_values = [1.0, 2.0]
    config = get_config()
    tmp_file = pwd() * tmp_dir * "tmp_params.jld2"

    #action
    KID.make_output_directories(tmp_dir)
    KID.save_data(u_values, u_names, config, file_name = tmp_file)
    tmp_data = KID.load_data(tmp_file)

    #test
    @test isdir(pwd() * tmp_dir)
    @test tmp_data.u_names == u_names
    @test tmp_data.u_bests == u_values
    @test tmp_data.config == config

    rm(pwd() * tmp_dir, recursive = true, force = true)
end
