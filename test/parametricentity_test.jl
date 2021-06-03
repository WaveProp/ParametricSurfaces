using Test
using ParametricSurfaces
using StaticArrays
using ParametricSurfaces: WavePropBase.HyperRectangle
using ParametricSurfaces: jacobian
using ParametricSurfaces: normal

@testset "curves" begin
    f   = (x) -> SVector(cos(x[1]),sin(x[1]))
    d   = HyperRectangle(0,2π)
    ent = ParametricEntity(f,d)
    s  =  SVector(rand())
    @test ent(s) == f(s)
    jac = jacobian(ent,s)
    @test jac[1] ≈ -sin(s[1]) && jac[2] ≈ cos(s[1])
    @test ent(s) ≈ normal(ent,s)
end

@testset "surfaces" begin
    f      = (x) -> SVector(x[1],x[2],sin(x[1]))
    d      = HyperRectangle((0,0),(1,1))
    ent    = ParametricEntity(f,d)
    s      =  SVector{2}(rand(2))
    @test ent(s) == f(s)
    @test jacobian(ent,s) ≈ [1 0; 0 1; cos(s[1]) 0]
end
