function gausslegendre(p,origin,width)
    nodes, weights = gausslegendre(p)
    @. nodes       = nodes*width/2  # shift to [-w/2,w/2]
    @. nodes       = nodes + width/2 + origin #shift to [a,b]
    @. weights     = weights * width/2
    return nodes, weights
end

function gausslegendre(p,hr::HyperRectangle{N}) where {N}
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

# if passed a single value of p, assume the same in all dimensions
gausslegendre(p::Integer,hr::HyperRectangle{N}) where {N} = gausslegendre(ones(Integer,N)*p,hr)
