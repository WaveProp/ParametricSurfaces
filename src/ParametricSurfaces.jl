module ParametricSurfaces

using StaticArrays
using LinearAlgebra
using ForwardDiff # for computing derivatives of parametric elements
using RecipesBase
using WavePropBase

import WavePropBase:
    AbstractEntity,
    AbstractElement,
    Domain,
    GenericMesh,
    ElementIterator,
    HyperRectangle,
    UniformCartesianMesh,
    ReferenceLine,
    ReferenceSquare,
    return_type,
    ambient_dimension,
    geometric_dimension,
    jacobian,
    new_tag,
    clear_entities!,
    global_add_entity!,
    center,
    domain,
    entities,
    assert_concrete_type,
    mesh

export
    #types
    ParametricEntity,
    HyperRectangle,
    #functions
    line

include("parametricentity.jl")
include("parametricelement.jl")
include("meshgen.jl")
include("simpleshapes.jl")

end # module
