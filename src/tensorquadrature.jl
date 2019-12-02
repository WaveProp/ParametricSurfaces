"""
Tensor quadrature for `M` dimensional space lifted to `N` dimensional space.

Contains also the normal vector at the quadrature nodes. The field `dims` contains the number of nodes in each dimension
. This can be used to recover the structured quadrature as
`reshape(quad,quad.dims...,:)`, which will yield `n` tensor quadratures of size `dims`
"""
struct TensorQuadrature{N,M,T}
    nodes::Vector{Point{N,T}}
    normals::Vector{Normal{N,T}}
    weights::Vector{T}
    dims::NTuple{M,Int}
end

Base.length(q::TensorQuadrature) = length(q.weights)

## Entity quadrature
function TensorQuadrature(p,surf::ParametricEntity{N,M,T},algo=gausslegendre) where {N,M,T}
    nelements = length(elements(surf))
    nnodes_per_element = prod(p)
    nnodes       = nnodes_per_element*nelements

    nodes      = Vector{Point{N,T}}(undef,nnodes)
    normals    = Vector{Normal{N,T}}(undef,nnodes)
    weights    = Vector{T}(undef,nnodes)

    n = 1
    for element in  elements(surf)
        _nodes, _weights = algo(p,element) #quadrature in reference element
        for (node,weight) in zip(_nodes,_weights)
            nodes[n] = surf(node)
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
            n+=1
        end
    end
    return TensorQuadrature{N,M,T}(nodes,normals,weights,Tuple(p))
end
# if passed a single value of p, assume the same in all dimensions
TensorQuadrature(p::Integer,surf::ParametricEntity{N,M},args...) where {N,M} = TensorQuadrature(ones(Integer,M)*p,surf,args...)


function TensorQuadrature(p,bdy::AbstractParametricBody{N,M,T},algo=gausslegendre) where {N,M,T}
    nodes_per_element = prod(p)
    nnodes = mapreduce(+,parts(bdy)) do part
        nelements = part |> elements |> length
        nodes_per_element*nelements
    end

    nodes      = Vector{Point{N,T}}(undef,nnodes)
    normals    = Vector{Normal{N,T}}(undef,nnodes)
    weights    = Vector{T}(undef,nnodes)

    n = 1
    for surf in parts(bdy)
        for element in  elements(surf)
            _nodes, _weights = algo(p,element) #quadrature in reference element
            for (node,weight) in zip(_nodes,_weights)
                nodes[n] = surf(node)
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
                n+=1
            end
        end
    end
    return TensorQuadrature{N,M,T}(nodes,normals,weights,Tuple(p))
end

################################################################################
@recipe function f(quad::TensorQuadrature{3,2})
    nodes = reshape(quad.nodes,quad.dims...,:)
    legend --> false
    grid   --> false
    aspect_ratio --> :equal
    seriestype := :surface
    # color  --> :blue
    linecolor --> :black
    # all patches
    for n =1:size(nodes,3)
        @series begin
            pts = nodes[:,:,n]
            x = [pt[1] for pt in pts]
            y = [pt[2] for pt in pts]
            z = [pt[3] for pt in pts]
            x,y,z
        end
    end
end
