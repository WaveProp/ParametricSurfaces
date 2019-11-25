"""
    ParametricEntity{Tf,Tx}

# Fields
- `f::Tf`
- `X::Tx`

Surface represented by image of `f(X)`, i.e. the image of `f`.
"""
struct ParametricEntity{F,Tx,Tq}
    f::F
    domain::Tx
end
(surf::ParametricEntity)(x) = surf.f(x)

jacobian(surf::ParametricEntity,s::AbstractArray) = jacobian(surf.f,s)
jacobian(surf::ParametricEntity,s::Point) = jacobian(surf.f,s)
jacobian(surf::ParametricEntity,s...)             = jacobian(surf.f,[s...])

function split(surf::ParametricEntity,axis,location)
    f = surf.f
    rec1, rec2 = split(surf.domain)
    return ParametricEntity(f,rec1), ParametricEntity(f,rec2)
end
