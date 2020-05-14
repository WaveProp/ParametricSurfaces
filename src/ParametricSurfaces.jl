module ParametricSurfaces

using  LinearAlgebra
using  RecipesBase
using  GeometryTypes: Point, Normal, HyperRectangle, Vec, vertices
using  GmshTools

import ForwardDiff: jacobian
import FastGaussQuadrature: gausschebyshev, gausshermite, gaussjacobi, gausslaguerre, gausslegendre, gausslobatto, gaussradau

export
    #types
    Circle,
    Ellipsis,
    Kite,
    Ellipsoid,
    Sphere,
    Bean,
    TensorQuadrature,
    GmshParametricEntity,
    GmshParametricBody,
    ParametricEntity,
    ParametricBody,
    #functions
    refine!,
    meshgen!,
    quadgen

include("hyperrectangle.jl")
include("parametricentity.jl")
include("parametricbody.jl")
include("tensorquadrature.jl")

end # module
