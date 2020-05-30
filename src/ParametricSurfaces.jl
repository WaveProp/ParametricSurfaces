module ParametricSurfaces

using  LinearAlgebra
using  RecipesBase
using  GeometryTypes: Point, Normal, HyperRectangle, Vec, vertices
using  GmshTools

import ForwardDiff

export
    #types
    Circle,
    Ellipsis,
    Kite,
    Ellipsoid,
    Sphere,
    Bean,
    Acorn,
    Cushion,
    TensorQuadrature,
    GmshParametricEntity,
    GmshParametricBody,
    ParametricEntity,
    ParametricBody,
    #functions
    refine!,
    meshgen!,
    quadgen

include("fejer.jl")
include("hyperrectangle.jl")
include("parametricentity.jl")
include("parametricbody.jl")
include("tensorquadrature.jl")

end # module
