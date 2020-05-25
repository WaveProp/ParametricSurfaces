abstract type AbstractEntity{N,M,T} end

"""
    ParametricEntity{N,M,T}

Represent an `M` dimensional surface `Y` embedded in `R^N` through the function
`f`; i.e. `Y` is given by `f : X -> Y`.
"""
struct ParametricEntity{N,M,T} <: AbstractEntity{N,M,T}
    parametrization::Function
    domain::HyperRectangle{M,T}
    elements::Vector{HyperRectangle{M,T}}
end

struct GmshParametricEntity{M} <: AbstractEntity{3,M,Float64}
    # dim=M
    tag::Int
    model::String
    domain::HyperRectangle{M,Float64}
    elements::Vector{HyperRectangle{M,Float64}}
end

getelements(ent::AbstractEntity) = ent.elements

" Dimension of Euclidean space where the entity is embedded"
ambient_dim(::ParametricEntity{N}) where {N} = N

"""
Geometrical dimension of entity; i.e. the number of parameters needed to describe it locally.

# line --> 1
# surface --> 2
# volume --> 3
"""
geo_dim(::ParametricEntity{N,M}) where {N,M} = M

function ParametricEntity(f,domain,els::Vector{<:HyperRectangle{M,T}}) where {M,T}
    pt = center(domain)
    N  = length(f(pt))
    ParametricEntity{N,M,T}(f,domain,els)
end

function GmshParametricEntity(dim::Int,tag::Int,model=gmsh.model.getCurrent())
    (umin,vmin),(umax,vmax) = gmsh.model.getParametrizationBounds(dim,tag)
    rec = HyperRectangle(umin,vmin,umax-umin,vmax-vmin)
    return GmshParametricEntity{dim}(tag,model,rec,[rec])
end
GmshParametricEntity(dim::Integer,tag::Integer,args...;kwargs...) = GmshParametricEntity(Int(dim),Int(tag),args...;kwargs...)

(par::ParametricEntity)(x) = par.parametrization(x)

function (par::GmshParametricEntity{N})(x) where {N}
    if N === 1
        return gmsh.model.getValue(N,par.tag,x)
    elseif N===2
        return gmsh.model.getValue(N,par.tag,[x[1],x[2]])
    else
        error("got N=$N, values must be 1 or 2")
    end
end

jacobian(psurf::ParametricEntity,s::AbstractArray) = jacobian(psurf.parametrization,s)
jacobian(psurf::ParametricEntity,s)                = jacobian(psurf.parametrization,[s...])
jacobian(psurf::ParametricEntity,s::Point)         = jacobian(psurf.parametrization,s)

function jacobian(psurf::GmshParametricEntity{N},s::Point) where {N}
    if N==1
        jac = gmsh.model.getDerivative(N,psurf.tag,s)
        return reshape(jac,3,N)
    elseif N==2
        jac = gmsh.model.getDerivative(N,psurf.tag,[s[1],s[2]])
        return reshape(jac,3,N)
    else
        error("got N=$N, values must be 1 or 2")
    end
end

function normal(ent::AbstractEntity,s)
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
function refine!(surf::AbstractEntity,ielem,axis)
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
function refine!(surf::AbstractEntity{N,M},ielem) where {N,M}
    if M == 1
        refine!(surf,ielem,1)
    elseif M == 2
        refine!(surf,ielem,1)
        n = length(getelements(surf))
        refine!(surf,ielem,2)
        refine!(surf,n,2)
    else
        @error "method not implemented"
    end
    return surf
end

#refine all elements in all directions
function refine!(surf::AbstractEntity)
    nel = length(getelements(surf))
    for i=1:nel
        refine!(surf,i)
    end
    return surf
end

function meshgen!(ent::AbstractEntity{N,M},h) where {N,M}
    domain   = ent.domain
    n        = ceil.(Int,(domain.widths)/h)
    elements = ent.elements
    resize!(elements,prod(n))
    if M == 1
        x = range(domain.origin[1],domain.origin[1]+domain.widths[1],length=n[1]+1)
        for i in 1:n[1]
            elements[i] = HyperRectangle(x[i],x[i+1]-x[i])
        end
    elseif M == 2
        x = range(domain.origin[1],domain.origin[1]+domain.widths[1],length=n[1]+1)
        y = range(domain.origin[2],domain.origin[2]+domain.widths[2],length=n[2]+1)
        cc = 1
        for i in 1:n[1]
            for j in 1:n[2]
                elements[cc] = HyperRectangle(x[i],y[j],x[i+1]-x[i],y[j+1]-y[j])
                cc += 1
            end
        end
    else
        error("type parameter $M must be 1 or 2")
    end
    return ent
end

################################################################################
@recipe function f(ent::ParametricEntity{2,1})
    els = ent.elements
    par = ent.parametrization
    legend --> false
    grid   --> false
    # aspect_ratio --> :equal
    for el in els
        @series begin
            pts     = [par(v) for v in vertices(el)]
            x       = [pt[1] for pt in pts]
            y       = [pt[2] for pt in pts]
            x,y
        end
    end
end

@recipe function f(ent::ParametricEntity{3,2},h=0.2)
    legend --> false
    grid   --> false
    # aspect_ratio --> :equal
    seriestype := :surface
    domain = ent.domain
    xmin,ymin = domain.origin
    xmax,ymax   = domain.origin + domain.widths
    xrange = xmin:h:xmax
    yrange = ymin:h:ymax
    pts    = [ent.parametrization((x,y)) for x in xrange, y in yrange]
    x      =  [pt[1] for pt in pts]
    y      =  [pt[2] for pt in pts]
    z      =  [pt[3] for pt in pts]
    x,y,z
end
