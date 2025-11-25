using SemiAnalyticWECs
using Test

@testset "DispersionEquations.jl" begin include("TestDispersionEquation.jl") end
@testset "EigenModes1D.jl" begin include("TestEigenModes1D.jl") end
@testset "GreensFunctions.jl" begin include("TestGreensFunctions.jl") end