struct TensorQuadrature{Ptype,Ntype,Wtype,M}
    nodes::Array{Ptype,M}
    normals::Array{Ntype,M}
    weights::Array{Wtype,M}
end

"""
    TensorQuadrature{N,T}(ent::ParametricEntity,algo)

Build a tensor product quadrature for the entity using the `algo` for doing a quadrature of `ent.domain`.

`algo` should have the signature `algo(::typeof(ent.domain))` and return `(::Vector{Vector{Ptype}},::Vector{Vector{Wtype}})`  representing the one-dimensional quadrature nodes.
"""
function TensorQuadrature{N,T}(surf::ParametricEntity,algo) where {N,T}
    domain = surf.domain
    _nodes, _weights = algo(domain)

    nodes      = similar(_nodes,Point{N,T})
    normals    = similar(_nodes,Normal{N,T})
    weights    = similar(_weights,T)

    for (n,(node,weight)) in enumerate(zip(_nodes,_weights))
        nodes[n] = surf(node)
        jac      = jacobian(surf,node)
        if N==2
            jac_det    = norm(jac)
            normals[n] = Normal{N,T}([jac[2],-jac[1]]./jac_det)
        elseif N==3
            tmp = cross(jac[:,1],jac[:,2])
            jac_det = norm(tmp)
            normals[n] = Normal{N,T}(tmp./jac_det)
        end
        weights[n] = jac_det*weight
    end
    TensorQuadrature{N,T,ndims(nodes)}(nodes,normals,weights)
end

function TensorQuadrature{N,T}(surf::ParametricEntity,order) where {N,T}
    domain = surf.domain
    _nodes, _weights = gausslegendre(order,domain)

    nodes      = similar(_nodes,Point{N,T})
    normals    = similar(_nodes,Normal{N,T})
    weights    = similar(_weights,T)

    for (n,(node,weight)) in enumerate(zip(_nodes,_weights))
        nodes[n] = surf(node)
        jac      = jacobian(surf,node)
        if N==2
            jac_det    = norm(jac)
            normals[n] = Normal{N,T}([jac[2],-jac[1]]./jac_det)
        elseif N==3
            tmp = cross(jac[:,1],jac[:,2])
            jac_det = norm(tmp)
            normals[n] = Normal{N,T}(tmp./jac_det)
        end
        weights[n] = jac_det*weight
    end
    TensorQuadrature{N,T,ndims(nodes)}(nodes,normals,weights)
end

function gausslegendre(p,origin,width)
    nodes, weights = gausslegendre(p)
    @. nodes       = nodes*width/2  # shift to [-w/2,w/2]
    @. nodes       = nodes + width/2 + origin #shift to [a,b]
    @. weights     = weights * width/2
    return nodes, weights
end

function gausslegendre(order,hr::HyperRectangle{N}) where {N}
    @assert length(order) == N
    _nodes   = Vector{Vector{Float64}}(undef,N)
    _weights = Vector{Vector{Float64}}(undef,N)
    for (n,p) in enumerate(order)
        _nodes[n],_weights[n] = gausslegendre(p,hr.origin[n],hr.widths[n])
    end
    nodes = Array{Point{N,Float64},N}(undef,order...)
    for (n,node) in enumerate(Iterators.product(_nodes...))
        nodes[n] = Point(node)
    end
    weights = Array{Float64,N}(undef,order...)
    for (n,weight) in enumerate(Iterators.product(_weights...))
        weights[n] = prod(weight)
    end
    return nodes,weights
end
