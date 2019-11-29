struct TensorQuadrature{N,T,Dtype}
    nodes::Vector{Point{N,T}}
    normals::Vector{Normal{N,T}}
    weights::Vector{T}
    dims::Dtype
end

"""
    TensorQuadrature{N,T}(ent::ParametricEntity,algo)

Build a tensor product quadrature for the entity using the `algo` for doing a quadrature of `ent.domain`.

`algo` should have the signature `algo(::typeof(ent.domain))` and return `(::Vector{Vector{Ptype}},::Vector{Vector{Wtype}})`  representing the one-dimensional quadrature nodes.
"""
function TensorQuadrature{P,N,W,D}(surf::ParametricEntity,p) where {P,N,W,D}
    domain = surf.domain
    _nodes, _weights = gausslegendre(p,domain)

    nodes      = Vector{P}(undef,length(_nodes))
    normals    = Vector{N}(undef,length(_nodes))
    weights    = Vector{W}(undef,length(_weights))

    for (n,(node,weight)) in enumerate(zip(_nodes,_weights))
        nodes[n] = surf(node)
        jac      = jacobian(surf,node)
        if length(P)==2
            jac_det    = norm(jac)
            normals[n] = N([jac[2],-jac[1]]./jac_det)
        elseif length(P)==3
            tmp = cross(jac[:,1],jac[:,2])
            jac_det = norm(tmp)
            normals[n] = N(tmp./jac_det)
        end
        weights[n] = jac_det*weight
    end
    TensorQuadrature{P,N,W,D}(nodes,normals,weights,p)
end

function TensorQuadrature(geo::ParametricBody,p)
    surf = geo.parts[1]
    domain = surf.domain
    _nodes, _weights = gausslegendre(p,domain)

    nodes      = Vector{P}(undef,length(_nodes))
    normals    = Vector{N}(undef,length(_nodes))
    weights    = Vector{W}(undef,length(_weights))

    for (n,(node,weight)) in enumerate(zip(_nodes,_weights))
        nodes[n] = surf(node)
        jac      = jacobian(surf,node)
        if length(P)==2
            jac_det    = norm(jac)
            normals[n] = N([jac[2],-jac[1]]./jac_det)
        elseif length(P)==3
            tmp = cross(jac[:,1],jac[:,2])
            jac_det = norm(tmp)
            normals[n] = N(tmp./jac_det)
        end
        weights[n] = jac_det*weight
    end
    TensorQuadrature{P,N,W,D}(nodes,normals,weights,p)
end

function TensorQuadrature(geo::ParametricBody,p)
    quads = []
    for surf in geo.parts
        quad = TensorQuadrature(surf,p)
        push!(quads,quad)
    end
    return quads
end

function gausslegendre(p,origin,width)
    nodes, weights = gausslegendre(p)
    @. nodes       = nodes*width/2  # shift to [-w/2,w/2]
    @. nodes       = nodes + width/2 + origin #shift to [a,b]
    @. weights     = weights * width/2
    return nodes, weights
end

function gausslegendre(p,hr::HyperRectangle{N}) where {N}
    @assert length(p) == N
    _nodes   = Vector{Vector{Float64}}(undef,N)
    _weights = Vector{Vector{Float64}}(undef,N)
    for (n,p) in enumerate(p)
        _nodes[n],_weights[n] = gausslegendre(p,hr.origin[n],hr.widths[n])
    end
    nodes = Array{Point{N,Float64},N}(undef,p...)
    for (n,node) in enumerate(Iterators.product(_nodes...))
        nodes[n] = Point(node)
    end
    weights = Array{Float64,N}(undef,p...)
    for (n,weight) in enumerate(Iterators.product(_weights...))
        weights[n] = prod(weight)
    end
    return nodes,weights
end
