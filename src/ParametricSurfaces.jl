module ParametricSurfaces

using  LinearAlgebra
using  RecipesBase
using  GeometryTypes: Point, Normal, HyperRectangle, Vec
using  GmshTools

import ForwardDiff: jacobian
import FastGaussQuadrature: gausslegendre

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
    refine!

include("hyperrectangle.jl")
include("parametricentity.jl")
include("parametricbody.jl")
include("tensorquadrature.jl")
include("gausslegendre.jl")

end # module
