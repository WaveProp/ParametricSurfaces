"""
    TensorQuadrature{N,M,T}
"""
struct TensorQuadrature{N,T}
    nodes::Vector{Point{N,T}}
    normals::Vector{Normal{N,T}}
    weights::Vector{T}
    elements::Vector{Vector{Int}}
end
Base.length(q::TensorQuadrature) = length(q.weights)
getnodes(q::TensorQuadrature)    = q.nodes
getnormals(q::TensorQuadrature)  = q.normals
getweights(q::TensorQuadrature)  = q.weights
getelements(q::TensorQuadrature) = q.elements
getnodes(q::TensorQuadrature,I)   = q.nodes[I]
getnormals(q::TensorQuadrature,I) = q.normals[I]
getweights(q::TensorQuadrature,I) = q.weights[I]
getelements(q::TensorQuadrature,I)= q.elements[I]

nodetype(q::TensorQuadrature)   = eltype(nodes(q))
normaltype(q::TensorQuadrature) = eltype(normals(q))
weighttype(q::TensorQuadrature) = eltype(weights(q))

# Base.permute!(quad::TensorQuadrature,perm::Vector{Int}) = map(x->permute!(x,perm),(quad.nodes,quad.normals,quad.weights))

## Entity quadrature
function TensorQuadrature(p,surf::AbstractEntity{N,M,T},algo=gausslegendre) where {N,M,T}
    nel      = length(getelements(surf))
    nnodes   = prod(p)*nel
    nodes    = Vector{Point{N,T}}(undef,nnodes)
    normals  = Vector{Normal{N,T}}(undef,nnodes)
    weights  = Vector{T}(undef,nnodes)
    elements = Vector{Vector{Int}}(undef,nel)

    n     = 0
    iel   = 0
    for element in  getelements(surf)
        iel += 1
        elements[iel] = []
        _nodes, _weights = _tensor_quad(p,element,algo) #quadrature in reference element
        for (node,weight) in zip(_nodes,_weights)
            n += 1
            nodes[n] = surf(node)
            push!(elements[iel],n) # add node index to current element list
            jac      = jacobian(surf,node)
            if N==2
                jac_det    = norm(jac)
                normals[n] = [jac[2],-jac[1]]./jac_det
            elseif N==3
                tmp = cross(jac[:,1],jac[:,2])
                jac_det = norm(tmp)
                normals[n] = tmp./jac_det
            end
            weights[n] = jac_det*weight
        end
    end
    return TensorQuadrature{N,T}(nodes,normals,weights,elements)
end


function TensorQuadrature(p,bdy::AbstractParametricBody{N,M,T},algo=gausslegendre) where {N,M,T}
    nelements = mapreduce(+,parts(bdy)) do part
        part |> getelements |> length
    end
    nnodes    = prod(p)*nelements
    nodes     = Vector{Point{N,T}}(undef,nnodes)
    normals   = Array{Normal{N,T}}(undef,nnodes)
    weights   = Vector{T}(undef,nnodes)
    elements  = Vector{Vector{Int}}(undef,nelements)

    n   = 0
    iel = 0
    for surf in parts(bdy)
        for patch in  getelements(surf)
            iel += 1
            elements[iel] = []
            _nodes, _weights = _tensor_quad(p,patch,algo) #quadrature in reference element
            for (node,weight) in zip(_nodes,_weights)
                n += 1
                nodes[n] = surf(node)
                push!(elements[iel],n) # add node index to current element list
                jac      = jacobian(surf,node)
                if N==2
                    jac_det    = norm(jac)
                    normals[n] = [jac[2],-jac[1]]./jac_det
                elseif N==3
                    tmp = cross(jac[:,1],jac[:,2])
                    jac_det = norm(tmp)
                    normals[n] = tmp./jac_det
                end
                weights[n] = jac_det*weight
            end
        end
    end
    return TensorQuadrature{N,T}(nodes,normals,weights,elements)
end

# if passed a single value of p, assume the same in all dimensions
TensorQuadrature(p::Integer,surf::AbstractEntity{N,M},args...) where {N,M} = TensorQuadrature(ntuple(i->p,M),surf,args...)
TensorQuadrature(p::Integer,surf::AbstractParametricBody{N,M},args...) where {N,M} = TensorQuadrature(ntuple(i->p,M),surf,args...)

################################################################################
@recipe function f(quad::TensorQuadrature{2})
    nodes = quad.nodes
    legend --> false
    grid   --> false
    aspect_ratio --> :equal
    # seriestype := :line
    # color  --> :blue
    # loop over elements
    n = 0
    for el in getelements(quad)
        n += 1
        @series begin
            linecolor --> n
            pts = getnodes(quad,el)
            x = [pt[1] for pt in pts]
            y = [pt[2] for pt in pts]
            x,y
        end
    end
end

@recipe function f(quad::TensorQuadrature{3})
    nodes = quad.nodes
    legend --> false
    grid   --> false
    # aspect_ratio --> :equal
    seriestype := :surface
    # color  --> :blue
    linecolor --> :black
    # loop over elements
    n = 0
    for el in getelements(quad)
        n += 1
        @series begin
            fillcolor --> n
            pts = getnodes(quad,el)
            x = [pt[1] for pt in pts]
            y = [pt[2] for pt in pts]
            z = [pt[3] for pt in pts]
            x,y,z
        end
    end
end
