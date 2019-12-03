using ParametricSurfaces
using Test

@testset "Tensor quadrature tests" begin
    let
        using ParametricSurfaces: HyperRectangle, gausslegendre
        line = HyperRectangle(-1.0,2.0)
        nodes,weights = gausslegendre(10,line)
        @test sum(weights) ≈ 2
        square = HyperRectangle(-1.0,-1.0,2.0,2.0)
        nodes,weights = gausslegendre((10,10),square)
        @test sum(weights) ≈ 4
        cube = HyperRectangle(-1.0,-1.0,-1.0,2.0,2.0,2.0)
        nodes,weights = gausslegendre((10,10,5),cube)
        @test sum(weights) ≈ 8
    end
end


@testset "Parametric entity tests" begin
    @testset "curves" begin
        let
            using ParametricSurfaces: HyperRectangle, ParametricEntity, jacobian, TensorQuadrature, refine!, normal
            f(x) = [cos(x[1]),sin(x[1])]
            domain = HyperRectangle(-1.0,2.0)
            ent  = ParametricEntity(f,[domain])
            s  =  rand()
            @test ent(s) == f(s)
            @test ent(s) ≈ normal(ent,s)
            @test jacobian(ent,[s]) ≈ [-sin(s),cos(s)]
            quad = TensorQuadrature(10,ent)
            @test length(quad) == 10
            @test sum(quad.weights) ≈ 2
            @test quad.normals ≈ quad.nodes
            refine!(ent)
            quad = TensorQuadrature(10,ent)
            @test length(quad) == 20
            @test sum(quad.weights) ≈ 2
            @test quad.normals ≈ quad.nodes
        end
    end
    @testset "surfaces" begin
        let
            using ParametricSurfaces: HyperRectangle, ParametricEntity, jacobian, refine!, Point
            f(x) = [x[1],x[2],0] #[-1,1]×[-1,1] square embedded in 3d
            domain = HyperRectangle(-1.0,-1.0,2.0,2.0)
            ent  = ParametricEntity(f,[domain])
            s  =  Point{2}(rand(2))
            @test ent(s) == f(s)
            @test jacobian(ent,s) ≈ [1 0; 0 1; 0 0]
            quad = TensorQuadrature((10,10),ent)
            @test length(quad) == 100
            @test quad.nodes == TensorQuadrature(10,ent).nodes
            @test sum(quad.weights) ≈ 4
            refine!(ent)
            quad = TensorQuadrature((10,10),ent)
            @test length(quad) == 4*100
            @test sum(quad.weights) ≈ 4
            quad = TensorQuadrature((5,10),ent)
            @test sum(quad.weights) ≈ 4
        end
    end
end

@testset "Parametric body tests" begin
    @testset "Circle tests" begin
        let
            using ParametricSurfaces: Circle, refine!, TensorQuadrature
            r    = π
            center = (1,0)
            geo  = Circle(;radius=π,center=center)
            quad = TensorQuadrature(10,geo)
            @test sum(quad.weights) ≈ 2π*r
            refine!(geo)
            quad = TensorQuadrature((10,),geo)
            @test sum(quad.weights) ≈ 2π*r
        end
    end
    @testset "Sphere tests" begin
        let
            using ParametricSurfaces: Sphere, refine!, Point, TensorQuadrature
            r   = 2
            geo = Sphere(radius=r)
            quad = TensorQuadrature((10,10),geo)
            @test length(quad) == 600
            @test !(sum(quad.weights) ≈ 4π*r^2)
            refine!(geo)
            quad = TensorQuadrature((10,10),geo)
            @test length(quad) == 4*600
            @test sum(quad.weights) ≈ 4π*r^2
        end
    end
    @testset "Ellipsoid tests" begin
        let
            using ParametricSurfaces: Ellipsoid, refine!, Point, TensorQuadrature
            r       = 2
            paxis   = (r,r,r)
            geo = Ellipsoid(paxis=paxis)
            quad = TensorQuadrature((10,10),geo)
            @test length(quad) == 600
            @test !(sum(quad.weights) ≈ 4π*r^2)
            refine!(geo)
            quad = TensorQuadrature((10,10),geo)
            @test length(quad) == 4*600
            @test sum(quad.weights) ≈ 4π*r^2
        end
    end
end
