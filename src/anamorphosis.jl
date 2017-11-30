module Anam

function loglin(t)

    function trans(x)
        if x < t
            return log(x)
        else
            return x/t - 1 + log(t)
        end
    end

    function invtrans(y)
        if y < log(t)
            return exp(y)
        else
            return (y + 1 - log(t)) * t
        end
    end

    return trans,invtrans
end

end
