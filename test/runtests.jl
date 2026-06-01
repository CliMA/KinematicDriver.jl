using SafeTestsets

@safetestset "Unit Tests" begin
    include("unit_tests/unit_test.jl")
end

@safetestset "Initial Profile Tests" begin
    include("initial_condition_tests/initial_profile_test.jl")
end

@safetestset "Aqua Tests" begin
    include("aqua.jl")
end
