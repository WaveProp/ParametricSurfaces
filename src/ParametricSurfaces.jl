module ParametricSurfaces

using  LinearAlgebra
using  RecipesBase
using  GeometryTypes: Point, Normal, HyperRectangle, Vec

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
    #functions
    refine!

include("hyperrectangle.jl")
include("parametricentity.jl")
include("parametricbody.jl")
include("tensorquadrature.jl")
include("gausslegendre.jl")

end # module
