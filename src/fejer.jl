function fejer1(n::Integer)
    theta = [(2j-1)*π/(2n) for j=1:n]
    x = -cos.(theta)
    w = zero(x)
    for j in 1:n
        tmp = 0.0
        for l in 1:ceil(n/2)
            tmp+= 1/(4*l^2-1) *cos(2*l*theta[j])
        end
        w[j] = 2/n * (1 - 2*tmp)
    end
    return x,w
end
