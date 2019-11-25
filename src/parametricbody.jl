struct ParametricBody{Tq}
    parts::Vector{ParametricEntity}
    quad::Vector{Tq}
end
ParametricBody(parts) = ParametricBody(parts,[])

function compute_quadrature!(bdy::ParametricBody,algo)
    for part in bdy.parts
        quad[n] = algo(part)
    end
end
