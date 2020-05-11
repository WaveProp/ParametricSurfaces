abstract type AbstractParametricBody{N,M,T} end

parts(bdy::AbstractParametricBody) = bdy.parts

refine!(bdy,args...;kwargs...) = map(part -> refine!(part,args...;kwargs...),parts(bdy))

# generic case, given simply by its parts
struct ParametricBody{N,M,T} <: AbstractParametricBody{N,M,T}
    parts::Vector{ParametricEntity{N,M,T}}
end

struct GmshParametricBody{M} <: AbstractParametricBody{3,M,Float64}
    parts::Vector{GmshParametricEntity{M}}
end

function GmshParametricBody(dim,tag,model=gmsh.model.getCurrent())
    dimtags = gmsh.model.getBoundary((dim,tag))
    body    = GmshParametricBody{2}([])
    for dimtag in dimtags
        patch = GmshParametricEntity(Int(dimtag[1]),Int(dimtag[2]),model)
        push!(body.parts,patch)
    end
    return body
end

struct Circle{T} <: AbstractParametricBody{2,1,T}
    center::Point{2,T}
    radius::T
    parts::Vector{ParametricEntity{2,1,T}}
end

function Circle{T}(;center=zeros(2),radius=1) where {T}
    f          = (s) -> radius.*[cospi(s[1]),sinpi(s[1])]
    domain     = HyperRectangle(-1.0,2.0)
    ent        = ParametricEntity(f,[domain])
    return Circle{T}(center,radius, [ent])
end
Circle(args...;kwargs...) = Circle{Float64}(args...;kwargs...)

struct Ellipsis{T} <: AbstractParametricBody{2,1,T}
    center::Point{2,T}
    paxis::Vec{2,T}
    parts::Vector{ParametricEntity{2,1,T}}
end
Ellipsis(args...;kwargs...) = Ellipsis{Float64}(args...;kwargs...)

function Ellipsis{T}(;center=zeros(2),paxis=ones(2)) where {T}
    f          = (s) -> paxis.*[cospi(s[1]),sinpi(s[1])]
    domain     = HyperRectangle(-1.0,2.0)
    surf       = ParametricEntity(f,[domain])
    return Ellipsis{T}(center,paxis,[surf])
end

struct Ellipsoid{T} <: AbstractParametricBody{3,2,T}
    center::Point{3,T}
    paxis::Vec{3,T}
    parts::Vector{ParametricEntity{3,2,T}}
end

function Ellipsoid{T}(;center=zeros(3),paxis=ones(3)) where {T}
    nparts = 6
    domain = HyperRectangle(-1.,-1.,2.,2.)
    parts  = Vector{ParametricEntity}(undef,nparts)
    for id=1:nparts
        param(x)     = _ellipsoid_parametrization(x[1],x[2],id,paxis,center)
        parts[id]    = ParametricEntity(param,[domain])
    end
    return Ellipsoid{T}(center,paxis,parts)
end
Ellipsoid(args...;kwargs...) = Ellipsoid{Float64}(args...;kwargs...)

struct Sphere{T} <: AbstractParametricBody{3,2,T}
    center::Point{3,T}
    radius::T
    parts::Vector{ParametricEntity{3,2,T}}
end

function Sphere{T}(;center=zeros(3),radius=1) where {T}
    nparts = 6
    domain = HyperRectangle(-1.,-1.,2.,2.)
    parts  = Vector{ParametricEntity}(undef,nparts)
    for id=1:nparts
        param(x)     = _sphere_parametrization(x[1],x[2],id,radius,center)
        parts[id]    = ParametricEntity(param,[domain])
    end
    return Sphere{T}(center,radius,parts)
end
Sphere(args...;kwargs...) = Sphere{Float64}(args...;kwargs...)

struct Bean{T} <: AbstractParametricBody{3,2,T}
    center::Point{3,T}
    paxis::Vec{3,T}
    parts::Vector{ParametricEntity{3,2,T}}
end
function Bean{T}(;center=zeros(3),paxis=ones(3)) where {T}
    nparts = 6
    domain = HyperRectangle(-1.,-1.,2.,2.)
    parts  = Vector{ParametricEntity}(undef,nparts)
    for id=1:nparts
        param(x)     = _bean_parametrization(x[1],x[2],id,paxis,center)
        parts[id]    = ParametricEntity(param,[domain])
    end
    return Bean{T}(center,paxis,parts)
end
Bean(args...;kwargs...) = Bean{Float64}(args...;kwargs...)

################################################################################
################################################################################
################################################################################
function _cube_parametrization(u,v,id,paxis,center)
    if id==1
        x = [1.,u,v]
    elseif id==2
        x = [-u,1.,v];
    elseif id==3
        x = [u,v,1.];
    elseif id==4
        x =[-1.,-u,v];
    elseif id==5
        x = [u,-1.,v];
    elseif id==6
        x = [-u,v,-1.];
    end
    return center .+ paxis.*x
end

function _sphere_parametrization(u,v,id,rad=1,center=zeros(3))
    if id==1
        x = [1.,u,v]
    elseif id==2
        x = [-u,1.,v];
    elseif id==3
        x = [u,v,1.];
    elseif id==4
        x =[-1.,-u,v];
    elseif id==5
        x = [u,-1.,v];
    elseif id==6
        x = [-u,v,-1.];
    end
    return center .+ rad.*x./sqrt(u^2+v^2+1)
end

function _ellipsoid_parametrization(u,v,id,paxis=ones(3),center=zeros(3))
    x = _sphere_parametrization(u,v,id)
    return x .* paxis .+ center
end

function _bean_parametrization(u,v,id,paxis=one(3),center=zeros(3))
    x = _sphere_parametrization(u,v,id)
    a = 0.8; b = 0.8; alpha1 = 0.3; alpha2 = 0.4; alpha3=0.1
    x[1] = a*sqrt(1.0-alpha3*cospi(x[3])).*x[1]
    x[2] =-alpha1*cospi(x[3])+b*sqrt(1.0-alpha2*cospi(x[3])).*x[2];
    x[3] = x[3];
    return x .* paxis .+ center
end

################################################################################
@recipe function f(bdy::AbstractParametricBody)
    for patch in bdy.parts
        @series begin
            patch
        end
    end
end
