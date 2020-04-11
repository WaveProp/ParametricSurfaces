"""
    ParametricEntity{N,M,T,F}

Represent an `M` dimensional surface `Y` embedded in `R^N` through the function
`f::F`; i.e. `Y` is given by `F : X -> Y`.
"""
struct ParametricEntity{N,M,T,F}
    parametrization::F
    elements::Vector{HyperRectangle{M,T}}
end
elements(ent::ParametricEntity) = ent.elements

" Dimension of Euclidean space where the entity is embedded"
ambient_dim(::ParametricEntity{N}) where {N} = N

"""
Geometrical dimension of entity; i.e. the number of parameters needed to describe it locally.

# line --> 1
# surface --> 2
# volume --> 3
"""
geo_dim(::ParametricEntity{N,M}) where {N,M} = M

function ParametricEntity(f,els::Vector{<:HyperRectangle{M,T}}) where {M,T}
    pt = center(first(els))
    N = length(f(pt))
    ParametricEntity{N,M,T,typeof(f)}(f,els)
end

(par::ParametricEntity)(x) = par.parametrization(x)

jacobian(psurf::ParametricEntity,s::AbstractArray) = jacobian(psurf.parametrization,s)
jacobian(psurf::ParametricEntity,s)                = jacobian(psurf.parametrization,[s...])
jacobian(psurf::ParametricEntity,s::Point)         = jacobian(psurf.parametrization,s)

function normal(ent::ParametricEntity,s)
    N        = ambient_dim(ent)
    jac      = jacobian(ent,s)
    if N==2
        jac_det    = norm(jac)
        return [jac[2],-jac[1]]./jac_det
    elseif N==3
        tmp = cross(jac[:,1],jac[:,2])
        jac_det = norm(tmp)
        return tmp./jac_det
    end
end

#split an element in single direction
function refine!(surf::ParametricEntity,ielem,axis)
    param     = surf.parametrization
    elem      = surf.elements[ielem]
    mid_point = elem.origin[axis]+elem.widths[axis]/2
    elem1, elem2    = split(elem, axis, mid_point)
    surf.elements[ielem] = elem1
    push!(surf.elements,elem2)
    return  surf
end

#refine an element in all directions (could be done better) TODO: improve this
#as it is very brittle at the moment since it assumes that the one dimensional
#refine! will ad its two new elements by (a) replacing the element ielem and
#pushing one elements to the back. Thus the "hacky" version below when M=2
function refine!(surf::ParametricEntity{N,M},ielem) where {N,M}
    if M == 1
        refine!(surf,ielem,1)
    elseif M == 2
         refine!(surf,ielem,1)
        n = length(elements(surf))
        refine!(surf,ielem,2)
        refine!(surf,n,2)
    else
        @error "method not implemented"
    end
    return surf
end

#refine all elements in all directions
function refine!(surf::ParametricEntity)
    nel = length(elements(surf))
    for i=1: nel
        refine!(surf,i)
    end
    return surf
end
