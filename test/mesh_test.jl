using Test
using ParametricSurfaces
using ParametricSurfaces: geometric_dimension, meshgen, Domain, boundary

@testset "Disk" begin
    disk = ParametricSurfaces.Disk() # abstract entity
    Γ = boundary(disk) |> Domain
    M = meshgen(Γ,(10,))
    # plot(M,Γ)
    @test geometric_dimension(disk) == 2
end

@testset "Ball" begin
    ball = ParametricSurfaces.Ball() # abstract entity
    Γ = boundary(ball) |> Domain
    M = meshgen(Γ,(1,1))
    # plot(M,Γ)
    @test geometric_dimension(ball) == 3
end
