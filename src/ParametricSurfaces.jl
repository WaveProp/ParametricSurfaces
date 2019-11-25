module ParametricSurfaces

using  GeometryTypes: Point, Normal, HyperRectangle
import ForwardDiff: jacobian
import  FastGaussQuadrature: gausslegendre
using  LinearAlgebra


include("parametricentity.jl")
include("parametricshapes.jl")
include("parametricbody.jl")
include("quadrature.jl")

end # module
