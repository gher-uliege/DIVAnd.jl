module Anam


function notransform()
    function trans(x; position = ())
        return x
    end

    function invtrans(y; position = ())
        return y
    end

    return trans,invtrans
end

"""
    trans,invtrans = loglin(t; epsilon = 0.)

Provide the following transform `log(x + epsilon)` (for x < t) and its inverse.
Beyond the threshold `t` (x ≥ t), the function is extended linearly in a
continous way.

`trans`,`invtrans` are scalar functions such that for any `x` (x > epsilon),
`x == invtrans(trans(x))`.

For any array `X`, we have: `X == invtrans.(trans.(X))`.
"""
function loglin(t; epsilon = 0.)

    function trans(x; position = ())
        if x < t
            return log(x + epsilon)
        else
            # derivative
            # D log(x + epsilon) = 1/(x+ϵ)
            # value at thereshold (x=t):  a = log(t+ϵ)
            # slope at thereshold (x=t):  b = 1/(t+ϵ)
            # linear extension
            # a + b (x-t)

            return log(t + epsilon) + (x-t) / (t+epsilon)
        end
    end

    function invtrans(y; position = ())
        if y < log(t + epsilon)
            # y = log(x + epsilon)
            # exp(y) - epsilon = x
            return exp(y) - epsilon
        else
            # y = log(t + epsilon) + (x-t) / (t+epsilon)
            # (y - log(t + epsilon)) * (t+epsilon)  = x-t

            return (y - log(t + epsilon)) * (t+epsilon) + t
        end
    end

    return trans,invtrans
end


"""
    trans,invtrans = logit(; min = 0., max = 1.)

Provide the logit transform and its inverse. Per default the logit transform
maps values within the interval from 0 and 1. This can be changed with the
`min` and `max` parameters. Note that trans(min) = -∞ and trans(max) = +∞.
The use safety-margin might be necessary.
"""
function logit(; min = 0., max = 1.)

    function trans(x; position = ())
        xs = (x-min) / (max-min)
        y = log(xs/(1-xs))
        return y
    end

    function invtrans(y; position = ())
        xs = 1/(1+exp(-y))
        x = min + xs * (max-min)
        return x
    end

    return trans,invtrans
end


end
