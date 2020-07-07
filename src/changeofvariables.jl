function _w(t,p)
    return @. 2π*_v(t,p)^p / ( _v(t,p)^p + _v(2π-t,p)^p)
end

function _v(t,p)
    return @. (1/p - 1/2)*((π-t)/π)^3 + 1/p*((t-π)/π) + 1/2
end

function smoothen(u,p,loc=:neither)
    if loc === :neither
        return u
    elseif loc === :both
        return -1 + 1/π*_w(π*(u+1),p)
    elseif loc === :left
        return -1 + 2/π*_w(π/2*(u+1),p)
    elseif loc === :right
        return -3 + 2/π*_w(π + π/2*(u+1),p)
    else
        @error "unknown option"
    end
end

function smoothen(ent::T,p,loc=:both) where {T<:ParametricEntity{2,1}}
    par_new = (u) -> ent.parametrization(smoothen(u,p,loc))
    return T(par_new,ent.domain,ent.elements)
end

function smoothen!(bdy::T,p,loc=:both) where {T<:AbstractParametricBody{2,1}}
    parts = getparts(bdy)
    for n in 1:length(parts)
        parts[n] = smoothen(parts[n],p,loc)
    end
    return bdy
end
