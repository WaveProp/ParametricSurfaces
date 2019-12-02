module ParametricSurfaces

export quadgen, circle

using  LinearAlgebra
using  RecipesBase
using  GeometryTypes: Point, Normal, HyperRectangle, Vec

import ForwardDiff: jacobian
import  FastGaussQuadrature: gausslegendre

include("hyperrectangle.jl")
include("parametricentity.jl")
include("parametricbody.jl")
include("tensorquadrature.jl")
include("gausslegendre.jl")
# include("parametricshapes.jl")
# include("geometry.jl")
# include("quadrature.jl")
# include("parametricentity.jl")
# include("parametricshapes.jl")
# include("parametricbody.jl")
# include("quadrature.jl")

end # module
