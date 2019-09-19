


function unstagger(x::AbstractRange)
    s = step(x)
    return s * (0:length(x)) .+ (first(x) - 0.5 * s)
end

function unstagger(x)
    return [
        x[1] - 0.5 * (x[2] - x[1]),
        (0.5 * (x[1:end-1] + x[2:end]))...,
        x[end] + 0.5 * (x[end] - x[end-1]),
    ]
end


function findin(x::AbstractRange, c)

    i = floor(Int, (c - first(x)) / step(x)) + 1
    if 1 <= i < length(x)
        return i
    else
        return -1
    end
end

function findin(x, c)
    for i = 1:length(x)-1
        if x[i] <= c < x[i+1]
            return i
        end
    end
    return -1
end


function binning(x::NTuple{1,T}, xv, v) where {T}
    # unstaggered coordinate
    x_unstagger = unstagger(x[1])
    sz = size(x[1])

    vb = zeros(eltype(v), sz)
    count = zeros(Int, sz)

    for j = 1:length(v)
        i = findin(x_unstagger, xv[1][j])
        if i != -1
            vb[i] = vb[i] + v[j]
            count[i] += 1
        end
    end

    return vb ./ count, count, vb
end
