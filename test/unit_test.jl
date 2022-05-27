"""
    A place to add unit tests (TODO)
"""

@testset "initial condition at the surface" begin

    init = KiD.init_condition(Float64, params, 0.0)

    @test init.θ == 297.9
    @test init.qv == 0.015
end
