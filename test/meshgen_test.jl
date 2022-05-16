using Test
using ParametricSurfaces
import WavePropBase as WPB

@testset "Disk" begin
    disk = ParametricSurfaces.Disk() # abstract entity
    Γ = WPB.boundary(disk) |> WPB.Domain
    M = ParametricSurfaces.meshgen(Γ,(10,))
    # plot(M,Γ)
    @test WPB.entities(Γ) == WPB.boundary(disk)
    @test WPB.geometric_dimension(disk) == 2
end

@testset "Ball" begin
    ball = ParametricSurfaces.Sphere() # abstract entity
    Γ = WPB.boundary(ball) |> WPB.Domain
    M = ParametricSurfaces.meshgen(Γ,(1,1))
    # plot(M,Γ)
    @test WPB.entities(Γ) == WPB.boundary(ball)
    @test WPB.geometric_dimension(ball) == 3
end
