######################## 2D ##############################################
function ellipsis(paxis=ones(2),center=ones(2))
    f(s)       = paxis.*[cospi(s[1]),sinpi(s[1])]
    domain     = HyperRectangle(-1.0,2.0)
    surf       = ParametricEntity(f,[domain])
    return Geometry{2}([surf])
end
circle(rad=1,center=ones(2)) = ellipsis(rad*ones(2),center)

function kite(rad=1,center=ones(2);q=20)
    f(s) = rad.*[cospi(s[1]) + 0.65*cospi(2*s[1]) - 0.65,
                1.5*sinpi(s[1])]
    domain = HyperRectangle(-1.0,2.0)
    surf   = ParametricEntity{2}(f,[domain])
    return Geometry{2}([surf])
end

######################## 3D ##############################################
function cube(paxis=ones(3),center=zeros(3))
    nparts = 6
    domain = HyperRectangle(-1.,-1.,2.,2.)
    parts = ParametricEntity{2,3,Float64}[]
    for id=1:nparts
        param(x) = __cube_parametrization(x[1],x[2],id,paxis,center)
        patch = ParametricEntity{3}(param,[domain])
        push!(parts,patch)
    end
    return Geometry{3}(parts)
end

function ellipsoid(paxis=ones(3),center=zeros(3))
    nparts = 6
    domain = HyperRectangle(-1.,-1.,2.,2.)
    parts = ParametricEntity{2,3,Float64}[]
    for id=1:nparts
        param(x) = __ellipsoid_parametrization(x[1],x[2],id,paxis,center)
        patch = ParametricEntity{3}(param,[domain])
        push!(parts,patch)
    end
    return Geometry{3}(parts)
end
sphere(rad=1,center=zeros(3)) = ellipsoid(rad*ones(3),center)

function bean(paxis=ones(3),center=zeros(3))
    nparts = 6
    domain = HyperRectangle(-1.,-1.,2.,2.)
    parts  = ParametricEntity{2,3,Float64}[]
    for id=1:nparts
        param(x) = __bean_parametrization(x[1],x[2],id,paxis,center)
        patch    = ParametricEntity{3}(param,[domain])
        push!(parts,patch)
    end
    return Geometry{3}(parts)
end
