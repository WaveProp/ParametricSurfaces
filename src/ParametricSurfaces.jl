module ParametricSurfaces

using  LinearAlgebra
using  RecipesBase
using  GeometryTypes: Point, Normal, HyperRectangle, Vec
using  GmshTools

import ForwardDiff: jacobian
import FastGaussQuadrature: gausschebyshev, gausshermite, gaussjacobi, gausslaguerre, gausslegendre, gausslobatto, gaussradau

export
    #types
    Circle,
    Ellipsis,
    Ellipsoid,
    Sphere,
    Bean,
    TensorQuadrature,
    GmshParametricEntity,
    GmshParametricBody,
    #functions
    refine!,
    gausschebyshev,
    gausshermite,
    gaussjacobi,
    gausslaguerre,
    gausslegendre,
    gausslobatto,
    gaussradau

include("hyperrectangle.jl")
include("parametricentity.jl")
include("parametricbody.jl")
include("tensorquadrature.jl")

end # module
