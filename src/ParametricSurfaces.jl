module ParametricSurfaces

using StaticArrays
using LinearAlgebra
using ForwardDiff # for computing derivatives of parametric elements

using WavePropBase
using WavePropBase.Utils
using WavePropBase.Geometry
using WavePropBase.Interpolation
using WavePropBase.Mesh

WavePropBase.@import_interface

export
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
