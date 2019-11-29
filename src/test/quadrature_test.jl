using ParametricSurfaces, GeometryTypes, Test
using ParametricSurfaces: circle, TensorQuadrature, gausslegendre, jacobian, kite, cube

# circle
r = 1.2
geo = circle(r)
quad = TensorQuadrature(geo,10)
@test sum(quad.weights) ≈ 2*π*r

# cube
geo  = cube()
quad = TensorQuadrature(geo.parts[1],(10,10))
quad = TensorQuadrature(geo,(10,10))
sum(quad.weights for quad in quad)

quad = TensorQuadrature(geo,10)


rec = HyperRectangle(-1.,-1,2,2)
nodes, weights = gausslegendre((10,10),rec)

tmp = Array{Point{2,Float64}}(undef,10,20)
for (n,node) in enumerate(Iterators.product(nodes...))
    tmp[n] = Point(node)
end
