struct ParametricBody
    parts::Vector{ParametricEntity}
end

function compute_quadrature!(bdy::ParametricBody,algo)
    for part in bdy.parts
        quad[n] = algo(part)
    end
end
