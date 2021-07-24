module ParametricSurfaces

using StaticArrays
using LinearAlgebra
using ForwardDiff # for computing derivatives of parametric elements
using RecipesBase

using WavePropBase
using WavePropBase.Utils
using WavePropBase.Geometry
using WavePropBase.Interpolation
using WavePropBase.Mesh

WavePropBase.@import_interface

export
    # re-exported from WavePropBase
    clear_entities!,
    ElementaryEntity,
    Domain,
    skeleton,
    internal_boundary,
    external_boundary,
    HyperRectangle,
    ReferenceLine,
    ReferenceSquare,
    #types
    ParametricEntity,
    ParametricElement,
    #functions
    line,
    meshgen

include("parametricentity.jl")
include("parametricelement.jl")
include("meshgen.jl")
include("simpleshapes.jl")

WavePropBase.@export_interface

end # module
