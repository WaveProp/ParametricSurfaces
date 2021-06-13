module ParametricSurfaces

using LinearAlgebra
using StaticArrays
using RecipesBase
using ForwardDiff

using WavePropBase

# interface methods
import WavePropBase:
    ambient_dimension,
    geometric_dimension,
    jacobian

include("parametricentity.jl")
include("simpleshapes.jl")
include("parametricelement.jl")
include("mesh.jl")

export
    # re-exported from WavePropBase
    ElementaryEntity,
    Domain,
    clear_entities!,
    HyperRectangle,
    skeleton,
    internal_boundary,
    external_boundary,
    jacobian,
    normal,
    # types
    ParametricEntity,
    # methods
    line,
    meshgen,
    flip_normal

end # module
