function get_linear_fit(x, y)
    X = zeros(length(x), 2)
    X[:,1] .= 1
    X[:,2] .= x
    return X \ y
end

function polyfit(x, y, n)
    X = zeros(length(x), n+1)
    X[:,1] .= 1
    for i = 2:(n+1)
        X[:,i] .= x .^ (i-1)
    end
    return X \ y
end

function eval_poly(x, coef)
    n = 0
    out = 0.0
    for c in coef
        out += c * x^n
        n += 1
    end
    return out
end

function residual(xvals, yvals, coef)
    res = 0.0
    tot = 0.0
    avgy = sum(yvals) / length(yvals)
    for (x,y) in zip(xvals, yvals)
        res = (y - eval_poly(x, coef))^2
        tot = (y - avgy)^2
    end

    return 1 - res/tot
end