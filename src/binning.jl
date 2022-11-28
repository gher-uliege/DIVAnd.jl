


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


function binning!(gridx::Tuple, x, v, vb, count)
    gridx_unstagger = DIVAnd.unstagger.(gridx)
    sz = length.(gridx)

    @assert size(vb) == sz
    @assert size(count) == sz

    nout = 0
    n = length(gridx)

    for j = 1:length(v)
        ind =
            CartesianIndex(ntuple(l -> DIVAnd.findin(gridx_unstagger[l], x[l][j]), Val(n)))

        if checkbounds(Bool, vb, ind)
            vb[ind] = vb[ind] + v[j]
            count[ind] += 1
        end
    end

    for ind = 1:length(vb)
        vb[ind] /= count[ind]
    end

    return vb, count
end

function binning(gridx::NTuple{N,AbstractVector}, x, v) where N
    sz = length.(gridx)
    vb = zeros(eltype(v), sz)
    count = zeros(Int, sz)

    return binning!(gridx, x, v, vb, count)
end


function binning(gridx::NTuple{N,Union{AbstractVector,DIVAnd.AbstractTimeSelector}}, x, v) where N
    @assert isa(gridx[N],DIVAnd.AbstractTimeSelector)

    gridx_ = gridx[1:N-1]
    TS  = gridx[N]

    sz = length.(gridx)
    vb = zeros(eltype(v), sz)
    mean_v = zeros(eltype(v), sz)
    count = zeros(Int, sz)

    for n = 1:length(gridx[N])
        sel = DIVAnd.select(TS,n,x[N])
        x_sel = ntuple(i -> x[i][sel],N-1)

        binning!(gridx_, x_sel, v[sel],selectdim(mean_v,N,n),selectdim(count,N,n))
    end

    return mean_v,count
end
