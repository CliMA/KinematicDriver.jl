using Test

import CLIMAParameters
import Kinematic1D

const CP = CLIMAParameters
const KiD = Kinematic1D

struct EarthParameterSet <: CP.AbstractEarthParameterSet end
const params = EarthParameterSet()

FT = Float64

@testset "initial condition at the surface" begin

    init = KiD.init_condition(FT, params, 0.0)

    @test init.θ == 297.9
    @test init.qv == 0.015
end
