struct ElementQuadrature{N,T}
    nodes::Array{T,2}
    normals::Array{T,2}
    weights::Vector{T}
end

"quadrature rule for integrating over a geometry. Assumes elements are all the same"
struct Quadrature{N,T}
    nodes::Vector{Point{N,T}}
    normals::Vector{Normal{N,T}}
    weights::Vector{T}
    nodes_per_element::Int
end

function quadgen(geo::Geometry{N,T},p) where {N,T}
    quad = Quadrature{N,T}([],[],[],p^(N-1))
    for  surf in geo.patches
        for el in surf.elements
            el_nodes, el_normals, el_weights = quadgen(surf,el,p)
            push!(quad.nodes,el_nodes...)
            push!(quad.weights,el_weights...)
            push!(quad.normals,el_normals...)
        end
    end
    return quad
end

function quadgen(surf::ParametricEntity{M,N,T}, domain::HyperRectangle, p) where {M,N,T}
    ref_nodes, ref_weights = quadgen(domain,p) # reference nodes and weights
    npts     = p^(N-1)
    nodes    = Vector{Point{N,T}}(undef,npts)
    normals  = Vector{Normal{N,T}}(undef,npts)
    weights  = Vector{T}(undef,npts)

    n = 1
    for (ref_node,ref_weight) in zip(ref_nodes,ref_weights)
        nodes[n]   = surf(ref_node)
        jac         = jacobian(surf,ref_node)
        if N==2
            jac_det = norm(jac)
            normals[n] = Normal{N,T}([jac[2],-jac[1]]./jac_det)
        elseif N==3
            tmp = cross(jac[:,1],jac[:,2])
            jac_det = norm(tmp)
            normals[n] = Normal{N,T}(tmp./jac_det)
        end
        weights[n] = jac_det*prod(ref_weight)
        n += 1
    end
    return nodes,normals,weights
end

# quadrature for (flat) hyper-rectangle
function quadgen(domain::HyperRectangle{N},p) where {N}
    N == 1 && return quad_interval(p,domain)
    N == 2 && return quad_square(p,domain)
end

function quad_1d(p,origin,width)
    nodes, weights = gausslegendre(p)
    #nodes, weights = newtoncotes(p)
    @. nodes       = nodes*width/2  # shift to [-w/2,w/2]
    @. nodes       = nodes + width/2 + origin #shift to [a,b]
    @. weights = weights * width/2
    return nodes, weights
end

function quad_interval(p,hr::HyperRectangle{1,T}) where {T}
    nodes, weights = quad_1d(p,hr.origin,hr.widths)
    return Vector{T}(nodes), Vector{T}(weights)
end

function quad_square(p,hr::HyperRectangle{2,T}) where {T}
    nodes_x, weights_x  = quad_1d(p,hr.origin[1],hr.widths[1])
    nodes_y, weights_y  = quad_1d(p,hr.origin[2],hr.widths[2])
    nodes    = vec([(nx,ny) for nx in nodes_x, ny in nodes_y])
    weights  = vec([ wx*wy for wx in weights_x, wy in weights_y])
    return nodes, weights
end

function Base.permute!(quad::Quadrature,perm::AbstractVector)
    permute!(quad.nodes,perm)
    permute!(quad.normals,perm)
    permute!(quad.weights,perm)
    return quad
end
