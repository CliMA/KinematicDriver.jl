"""
    A place to add unit tests (TODO)
"""

@testset "initial condition at the surface" begin

    init = KiD.init_condition(Float64, params, 0.0)

    @test init.θ_std == 297.9
end
