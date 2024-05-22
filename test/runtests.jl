using SafeTestsets

@safetestset "Unit Tests" begin
    include("unit_tests/unit_test.jl")
end

@safeTestset "Optimization Tests" begin
    include("opt_tests/opt_tests.jl")
end

@safetestset "Initial Profile Tests" begin
    include("initial_condition_tests/initial_profile_test.jl")
end

@safetestset "Box Test" begin
    include("experiments/box_driver/run_box_simulation.jl")
end

@safetestset "K1D Driver Test" begin
    include("experiments/KiD_driver/KiD_driver.jl")
end

@safetestset "Collision Sedimentation Tests" begin
    include("experiments/KiD_col_sed_driver/run_KiD_col_sed_simulation.jl")
end

@safetestset "K2D Simulation Tests" begin
    include("experiments/Ki2D_driver/run_kinematic2d_simulations.jl")
end

@safetestset "Aqua Tests" begin
    include("aqua.jl")
end
