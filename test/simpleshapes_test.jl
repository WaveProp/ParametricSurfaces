using Test
using ParametricSurfaces
using ParametricSurfaces: geometric_dimension

@testset "Disk" begin
    disk = ParametricSurfaces.Disk() # abstract entity
    @test geometric_dimension(disk) == 2
end

@testset "Ball" begin
    ball = ParametricSurfaces.Ball() # abstract entity
    @test geometric_dimension(ball) == 3
end
