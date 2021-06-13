using Test
using ParametricSurfaces
using StaticArrays
using ParametricSurfaces: ParametricElement

@testset "Curves" begin
    d = HyperRectangle(-1,1)
    el = ParametricElement(d) do u
        s = u[1]
        SVector(s,s^2)
    end
    @test el(0)   == SVector(-1,1)
    @test el(1)   == SVector(1,1)
    @test el(0.5) == SVector(0,0)
end
