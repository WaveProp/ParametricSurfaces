# Some simple shapes coded as AbstractEntities for convenience
struct Disk <: AbstractEntity
    # dim = 2
    tag::Int
    center::SVector{2,Float64}
    paxis::SVector{2,Float64}
    boundary::Vector{ParametricEntity}
    function Disk(t,c,axis,bnd)
        ent = new(t,c,axis,bnd)
        global_add_entity!(ent)
        return ent
    end
end

boundary(d::Disk) = d.boundary

function Disk(;center=(0, 0),paxis=(1,1))
    f          = (s) -> paxis.*SVector(cospi(s[1]),sinpi(s[1]))
    domain     = HyperRectangle(-1.0,1.0)
    ent        = ParametricEntity(f, domain)
    tag        = new_tag(2) # generate a unique tag for entities of dimension 2
    return Disk(tag,center,paxis,[ent])
end

geometric_dimension(::Disk) = 2
tag(d::Disk)                = d.tag

struct Ball <: AbstractEntity
    # dim = 3
    tag::Int
    center::SVector{3,Float64}
    radius::Float64
    boundary::Vector{ParametricEntity}
    function Ball(t,c,r,bnd)
        ent = new(t,c,r,bnd)
        global_add_entity!(ent)
        return ent
    end
end

boundary(b::Ball) = b.boundary
geometric_dimension(::Ball) = 3
ambient_dimension(::Ball)   = 3
tag(d::Ball)                = d.tag

function Ball(;center=(0,0,0),radius=1)
    nparts = 6
    domain = HyperRectangle((-1.,-1.),(1.,1))
    parts  = Vector{ParametricEntity}(undef,nparts)
    for id=1:nparts
        param     = (x) -> _sphere_parametrization(x[1],x[2],id,radius,center)
        parts[id] = ParametricEntity(param,domain)
    end
    tag = new_tag(3)
    return Ball(tag,center,radius,parts)
end

# parametrize the sphere (surface) with 6 patches
function _sphere_parametrization(u,v,id,rad=1,center=zeros(3))
    # parametrization of 6 patches
    if id==1
        x = SVector(1.,u,v)
    elseif id==2
        x = SVector(-u,1.,v)
    elseif id==3
        x = SVector(u,v,1.)
    elseif id==4
        x = SVector(-1.,-u,v)
    elseif id==5
        x = SVector(u,-1.,v)
    elseif id==6
        x = SVector(-u,v,-1.)
    end
    return center .+ rad.*x./sqrt(u^2+v^2+1)
end
