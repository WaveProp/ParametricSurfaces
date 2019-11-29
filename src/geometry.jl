struct Geometry{N,T}
    patches::Vector{ParametricEntity}
end

Geometry{N}(args...) where {N} = Geometry{N,Float64}(args...)

refine!(geo::Geometry) = for patch in geo.patches; refine!(patch); end
