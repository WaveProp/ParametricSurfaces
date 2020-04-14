origin(rec::HyperRectangle) = rec.origin
widths(rec::HyperRectangle) = rec.widths
center(rec::HyperRectangle) = origin(rec) .+ widths(rec)

function _tensor_quad(p,origin,width,algo)
    nodes, weights = algo(p)
    @. nodes       = nodes*width/2  # shift to [-w/2,w/2]
    @. nodes       = nodes + width/2 + origin #shift to [a,b]
    @. weights     = weights * width/2
    return nodes, weights
end

function _tensor_quad(p,hr::HyperRectangle{N},algo) where {N}
    _nodes   = Vector{Vector{Float64}}(undef,N)
    _weights = Vector{Vector{Float64}}(undef,N)
    for (n,p) in enumerate(p)
        _nodes[n],_weights[n] = _tensor_quad(p,hr.origin[n],hr.widths[n],algo)
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

# if passed a single value of p, assume the same in all dimensions
_tensor_quad(p::Integer,hr::HyperRectangle{N}) where {N} = _tensor_quad(ones(Integer,N)*p,hr)
