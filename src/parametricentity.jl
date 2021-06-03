"""
    ParametricEntity{D,F} <: AbstractEntity

A geometric entity with an explicit `parametrization`. The entity is the image
of `domain` under the parametrization.
"""
struct ParametricEntity <: AbstractEntity
    dim::UInt8
    tag::Int
    parametrization
    domain
    function ParametricEntity(dim, tag, f, d)
        ent = new(dim, tag, f, d)
        # every entity gets added to a global variable ENTITIES so that we can
        # ensure the (d,t) pair is a UUID for an entity, and to easily retrieve
        # different entities.
        global_add_entity!(ent)
        return ent
    end
end

domain(p::ParametricEntity)              = p.domain
parametrization(p::ParametricEntity)     = p.parametrization
geometric_dimension(p::ParametricEntity) = p.dim
tag(p::ParametricEntity) = p.tag

# TODO: `boundary` method missing. One can always extract the boundary of the
# parametric entity by "manually" taking the boundary of domain and
# reparametrizing (top-down approach). Postponing the implementation of this
# untile I know how this could be useful.

function ParametricEntity(f,dom)
    d = geometric_dimension(dom)
    t = new_tag(d) # automatically generate a new (valid) tag
    ParametricEntity(d, t, f, dom)
end

function return_type(p::ParametricEntity)
    # NOTE: this relies on the brittle promote_op
    d = domain(p)
    x = center(d)
    f = parametrization(p)
    T = Base.promote_op(f, typeof(x))
    return T
end

ambient_dimension(p::ParametricEntity) = length(return_type(p))

function (par::ParametricEntity)(x)
    @assert x in domain(par)
    par.parametrization(x)
end

jacobian(psurf::ParametricEntity,s::SVector)  = ForwardDiff.jacobian(psurf.parametrization, s::SVector)
jacobian(psurf::ParametricEntity,s)           = jacobian(psurf, SVector(s))

function normal(ent::ParametricEntity, u)
    dim = geometric_dimension(ent)
    t   = tag(ent)
    @assert length(u) == dim
    s    = sign(t) # flip the normal if negative sign
    N    = ambient_dimension(ent)
    M    = geometric_dimension(ent)
    msg  = "computing the normal vector requires the entity to be of co-dimension one."
    @assert N - M == 1 msg
    if M == 1 # a line in 2d
        t = jacobian(ent, u)
        ν = SVector(t[2], -t[1])
        return s * ν / norm(ν)
    elseif M == 2 # a surface in 3d
        j  = jacobian(ent, u)
        t₁ = j[:,1]
        t₂ = j[:,2]
        ν  = cross(t₁, t₂)
        return s * ν / norm(ν)
    else
        notimplemented()
    end
end
