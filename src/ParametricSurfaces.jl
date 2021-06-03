module ParametricSurfaces

using LinearAlgebra
using StaticArrays
using RecipesBase
using ForwardDiff

using WavePropBase

@interface geometric_dimension
@interface ambient_dimension

include("parametricentity.jl")
include("simpleshapes.jl")
include("parametricelement.jl")
include("mesh.jl")

export
    # types
    ParametricEntity

end # module
