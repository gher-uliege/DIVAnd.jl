using Random: seed!

"""
    meanx,stdx = stats(sumx,sumx2,N)

Computes the mean `meanx` and the standard deviation `stdx` from the sum
(`sumx`) and the sum of squares (`sumx2`) from `N` numbers.
"""
function stats(sumx, sumx2, N)
    # var(x) = 1/(N-1)  Σᵢ (xᵢ - mean(x))²
    # var(x) = 1/(N-1)  Σᵢ [ xᵢ² - 2 xᵢ mean(x) + mean(x)² ]
    # var(x) = 1/(N-1) ( Σᵢ xᵢ²  -  Σᵢ 2 xᵢ mean(x) + N mean(x)² )
    # var(x) = 1/(N-1) ( Σᵢ xᵢ²  -  2 N mean(x)² + N mean(x)² )
    # var(x) = 1/(N-1) ( Σᵢ xᵢ²  -  N mean(x)² )

    meanx = sumx / N
    stdx = sumx2 - N * meanx .^ 2

    # de-bias std
    stdx = stdx / (N - 1)

    stdx = if stdx < 0
        zero(stdx)
    else
        sqrt(stdx)
    end

    return meanx, stdx
end

"""
    meanx,meany,stdx,stdy,covar,corr = stats(sumx,sumx2,sumy,sumy2,sumxy,N)

Computes the mean `meanx` and the standard deviation `stdx` from the sum
(`sumx`) and the sum of squares (`sumx2`) from `N` numbers and
similarly for the variable `y`. The function computes also the
Pearson correlation `corr` and covariance `covar` between `x` and `y`.
"""
function stats(sumx, sumx2, sumy, sumy2, sumxy, N)
    # (N-1) * covar(x,y) = Σᵢ (xᵢ - mean(x)) (yᵢ - mean(y))
    # (N-1) * covar(x,y) = Σᵢ (xᵢ yᵢ - xᵢ mean(y) - mean(x) yᵢ + mean(x) mean(y))
    # (N-1) * covar(x,y) = Σᵢ xᵢ yᵢ - N mean(x) mean(y) - N mean(x) mean(y) + N mean(x) mean(y)
    # (N-1) * covar(x,y) = Σᵢ xᵢ yᵢ - N mean(x) mean(y)

    meanx, stdx = stats(sumx, sumx2, N)
    meany, stdy = stats(sumy, sumy2, N)

    covar = sumxy - sumx * sumy / N

    covar = covar / (N - 1)

    corr = covar / (stdx * stdy)

    return meanx, meany, stdx, stdy, covar, corr
end


function fitchoose(x)
    # random index
    i = rand(1:length(x[1]))
    j = rand(1:length(x[1]))
    return (i, j)
end

"""
    distx,covar,corr,varx,count = empiriccovar(x,v,distbin,mincount;
                              maxpoints = 10000,
                              distfun = (xi,xj) -> sqrt(sum(abs2,xi-xj)))

Compute the covariance, correlation and variance of a cloud of data points with
the value `v` (a vector) and the location `x` (a tuple of vectors) grouped by
distance. Random pairs are choosen and grouped by their distance
(computed by `distfun`) in bins defined by `distbin`. The function try to fill
at least `mincount` of data points in each bin but always stop after
considering `maxpoints` pairs.
"""
function empiriccovar(
    x,
    v,
    distbin,
    mincount;
    maxpoints = 10000,
    choose = fitchoose,
    distfun = (xi, xj) -> sqrt(sum(abs2, xi - xj)),
)

    @assert all(length.(x) .== length(v))


    pmax = length(distbin) - 1
    sumvivj = zeros(pmax)
    sumvi = zeros(pmax)
    sumvj = zeros(pmax)
    sumvi2 = zeros(pmax)
    sumvj2 = zeros(pmax)
    count = zeros(pmax)

    corr = zeros(pmax)
    covar = zeros(pmax)
    varx = zeros(pmax)
    count = zeros(pmax)

    distx = zeros(pmax)
    for i = 1:pmax
        distx[i] = (distbin[i] + distbin[i+1]) / 2
    end

    # coordinates of considered points
    xi = zeros(eltype(x[1]), length(x))
    xj = zeros(eltype(x[1]), length(x))

    for l = 1:maxpoints
        # random index
        i, j = choose(x)

        for k = 1:length(x)
            xi[k] = x[k][i]
            xj[k] = x[k][j]
        end

        if isnan(v[i]) || isnan(v[j])
            # one point is masked
            continue
        end

        distance = distfun(xi, xj)

        if distance >= distbin[end]
            # distance too large
            continue
        end

        p = findlast(distance .>= distbin)

        #if count[p] >= mincount
        #    # already enought points
        #    continue
        #end

        #distbin(p) <= dist && dist < distbin(p+1)

        sumvivj[p] += v[i] * v[j]
        sumvi[p] += v[i]
        sumvj[p] += v[j]
        sumvi2[p] += v[i]^2
        sumvj2[p] += v[j]^2
        count[p] += 1

        #if mod(l,1000) == 0
        #    @show count
        #end

        if all(count .>= mincount)
            break
        end
    end

    for i = 1:pmax
        meanvi, meanvj, stdvi, stdvj, covar[i], corr[i] =
            stats(sumvi[i], sumvi2[i], sumvj[i], sumvj2[i], sumvivj[i], count[i])
        varx[i] = ((stdvi^2 + stdvj^2) * count[i]) / (2 * count[i] - 1)

        #meanx,stdx = stats(sumvi[i] + sumvj[i],sumvi2[i] + sumvj2[i],2*count[i])
        #@show i, stdx, sqrt(varx[i])
        #@show i,distbin[i],distbin[i+1],covar[i],corr[i],sumvivj[i],count[i]
    end


    return distx, covar, corr, varx, count
end



# mean over dimensions 2 ignoring NaNs
function nm(covar::Array{T,2}) where {T}
    m = zeros(T, size(covar, 1))
    s = zeros(T, size(covar, 1))
    count = zeros(Int, size(covar, 1))
    for j = 1:size(covar, 2)
        for i = 1:size(covar, 1)
            if !isnan(covar[i, j])
                m[i] = m[i] + covar[i, j]
                s[i] = s[i] + covar[i, j]^2
                count[i] += 1
            end
        end
    end

    for i = 1:size(covar, 1)
        m[i], s[i] = stats(m[i], s[i], count[i])
    end

    return m, s
end

function empiriccovarmean(
    x,
    v::Vector{T},
    distbin,
    mincount;
    maxpoints = 10000,
    nmean = 10,
    choose = fitchoose,
    distfun = (xi, xj) -> sqrt(sum(abs2, xi - xj)),
)::NTuple{6,Vector{T}} where {T}


    sz = (length(distbin) - 1, nmean)
    covar = zeros(sz)
    corr = zeros(sz)
    varx = zeros(sz)
    count = zeros(sz)

    # make sure that the variable distx is visible
    # outside the for loop
    # https://web.archive.org/save/https://docs.julialang.org/en/release-0.6/manual/variables-and-scoping/#Soft-Local-Scope-1
    # https://web.archive.org/web/20180217105449/https://stackoverflow.com/questions/22798305/function-variable-does-not-live-outside-a-for-loop

    # necessary for type inference of distx

    distx, covar[:, 1], corr[:, 1], varx[:, 1], count[:, 1] = empiriccovar(
        x,
        v,
        distbin,
        mincount;
        maxpoints = maxpoints,
        choose = choose,
        distfun = distfun,
    )


    for k = 2:nmean
        distx, covar[:, k], corr[:, k], varx[:, k], count[:, k] = empiriccovar(
            x,
            v,
            distbin,
            mincount;
            maxpoints = maxpoints,
            choose = choose,
            distfun = distfun,
        )
    end

    meancovar, stdcovar = nm(covar)
    meancorr, stdcorr = nm(corr)
    meanvarx, stdvarx = nm(varx)
    meancount, stdcount = nm(count)

    #return distx,nm(covar),nm(corr), nm(varx), nm(count)
    return distx::Vector{T}, meancovar, meancorr, meanvarx, meancount, stdcovar
    #    stdcovar,meancount
end


function distfun_euclid(x0, x1)
    dist = zero(x0[1])
    for i = 1:length(x0)
        dist += (x0[i] - x1[i])^2
    end
    return dist = sqrt(dist)
end


distfun_m(x0, x1) = EarthRadius * distance(x0[2], x0[1], x1[2], x1[1]) * pi / 180


mutable struct AllCoupels
    n::Int
end

function Base.iterate(iter::AllCoupels, state = (1, 1))
    i, j = state
    if (i == iter.n - 1) && (j == iter.n)
        return nothing
    end

    if j < iter.n
        nextstate = (i, j + 1)
    else
        nextstate = (i + 1, i + 2)
    end

    return (nextstate, nextstate)
end

Random.seed!(iter::AllCoupels,iseed) = nothing

mutable struct RandomCoupels{TRNG <: AbstractRNG}
    n::Int
    count::Int
    rng::TRNG
end

Base.length(iter::RandomCoupels) = iter.count
Random.seed!(iter::RandomCoupels,iseed) = Random.seed!(iter.rng,iseed)

function Base.iterate(iter::RandomCoupels, state = (0, copy(iter.rng)))
    count, rng = state

    if count == iter.count
        return nothing
    end

    # pick two random points
    j = rand(rng, 1:iter.n)
    i = j
    while (i == j)
        i = rand(rng, 1:iter.n)
    end

    return ((i, j), (count + 1, rng))
end

mutable struct VertRandomCoupels{TRNG <: AbstractRNG,Tz,Tx,Ts}
    zlevel::Tz # depth in meters
    zindex::Vector{Int}
    x::NTuple{3,Vector{Tx}}
    searchxy::Ts # in meters
    maxntries::Int
    count::Int
    rng::TRNG
end

Random.seed!(iter::VertRandomCoupels,iseed) = Random.seed!(iter.rng,iseed)

function _next(iter::VertRandomCoupels, state)
    count, rng = state

    # pick two random points
    j = -1
    jindex = -1

    for ntries = 1:iter.maxntries
        j = iter.zindex[rand(iter.rng,1:length(iter.zindex))]

        # # for few data, this is faster
        # if length(zindex) < 10000
        #     mask = falses(size(x[1]))
        #     for k = 1:length(x[1])
        #         mask[k] = distfun([x[1][k],x[2][k]],[x[1][j],x[2][j]]) < searchxy
        #     end
        #     jindex = pickone(mask)
        #     return j,jindex
        # end

        jindex = -1

        for ntries2 = 1:iter.maxntries
            k = rand(iter.rng,1:length(iter.x[1]))
            if (
                (
                    distfun_m([iter.x[1][k], iter.x[2][k]], [iter.x[1][j], iter.x[2][j]]) <
                    iter.searchxy
                ) && (j != k)
            )
                #            if (distfun_m([iter.x[1][k],iter.x[2][k]],[iter.x[1][j],iter.x[2][j]]) < iter.searchxy)
                jindex = k
                break
            end
        end

        if jindex != -1
            break
        end
    end

    if (j == -1) || (jindex == -1)
        error("fail to find enought pairs at z = $(iter.zlevel)")
    end
    #@show iter.zlevel,iter.x[3][j],iter.x[3][jindex]
    return ((j, jindex), (count + 1, rng))
end

Base.iterate(iter::VertRandomCoupels, state = (0, copy(iter.rng))) =
    (state[1] == iter.count ? nothing : _next(iter, state))



function fitlen(x::Tuple, d, nsamp; kwargs...)
    weight = ones(size(d))
    return fitlen(x, d, weight, nsamp; kwargs...)
end

"""
  varbak, RL, dbinfo =
     fitlen(x::Tuple,d,weight,nsamp; distfun = distfun_euclid, kwargs...)

This function used to be called lfit in fitlsn.f

"""
function fitlen(x::Tuple, d, weight, nsamp;
                rng = Random.GLOBAL_RNG, kwargs...)
    # number of samples
    n = length(d)
    iseed = n
    Random.seed!(rng,iseed)

    iter = if (nsamp == 0)
        AllCoupels(n)
    else
        @debug "will generate random couples"
        if (nsamp > n)
            @warn "Strange to ask for more samples than available from data; will proceed"
        end

        RandomCoupels(n, (nsamp * (nsamp - 1)) ÷ 2, rng)
    end

    #    if (n > 10000) && (nsamp != 0)
    #        @warn "Be patient big data set: $n"
    #    end

    return fitlen(x::Tuple, d, weight, nsamp, iter; kwargs...)
end


function fitlen(
    x::Tuple,
    d,
    weight,
    nsamp,
    iter;
    distfun = distfun_euclid,
    iseed = length(d),
    kwargs...,
)
    if length(d) == 0
        @warn "no data is provided to fitlen"
        dbinfo = Dict{Symbol,Any}(
        :covar => [],
        :fitcovar => [],
        :distx => [],
        :rqual => 0.0,
        :range => 1:0,
        :covarweight => [],
        :distx => [],
        :sn => 0.0,
        :meandist => 0.0,
    )
        return NaN, NaN, dbinfo
    end

    # number of dimensions
    ndims = length(x)

    # number of samples
    n = length(d)

    # per default operate on all data
    nop = n

    rqual = 0.0
    maxdist = 0.0
    meandist = 0.0
    dist = 0.0
    rjjj = 0.0


    # compute mean and variance using the weights
    datamean = 0.0
    datavar = 0.0
    rn = 0.0

    for i = 1:n
        datamean = datamean + d[i] * weight[i]
        datavar = datavar + d[i] * d[i] * weight[i]
        rn = rn + weight[i]
    end

    datamean = datamean / rn
    variance = datavar / rn - datamean^2
    @debug "Number of data points: $n"
    @debug "data mean: $datamean"
    @debug "data variance: $variance"

    @debug "Now calculating distance distribution"

    x0 = zeros(ndims)
    x1 = zeros(ndims)

    Random.seed!(iter,iseed)

    for (i, j) in iter
        # compute the distance
        for l = 1:ndims
            x0[l] = x[l][i]
            x1[l] = x[l][j]
        end
        dist = distfun(x0, x1)

        meandist = meandist + dist
        if (dist > maxdist)
            maxdist = dist
        end
    end

    if (nsamp == 0)
        rjjj = rn * (rn - 1.0) * 0.5
    else
        rjjj = nsamp * (nsamp - 1) / 2
    end

    @debug "Number of data couples considered: $rjjj"
    meandist = meandist / rjjj

    @debug "maximum distance between points: $maxdist"

    @debug "Mean distance between points: $meandist"

    rnbins = if (nsamp == 0)
        min(80.0, rn^2 / maxdist * meandist / 20.0)
    else
        min(80.0, nsamp^2 / maxdist * meandist / 20.0)
    end

    @debug "Number of probable active bins: $rnbins"

    ddist = meandist / rnbins
    nbmax=1
    worktmp=maxdist / ddist
    if !isnan(worktmp)&&isfinite(worktmp)
             nbmax = floor(Int, worktmp + 1)
    end


    #nbmax = floor(Int, maxdist / ddist + 1)
    @debug "distance for binning: $ddist"
    @debug "maximum number of bins: $nbmax"

    if (nsamp == 0)
        @debug "Average number of pairs in each bin: $(rn*rn/nbmax/2)"
    else
        @debug "Average number of pairs in each bin: $(nsamp*nsamp/nbmax/2)"
    end

    # d_i' are anomalies (d_i - datamean)

    # sum of d_i' * d_j' * w_i  * w_j
    covar = zeros(nbmax)

    # sum of (d_i' * d_j')² * w_i  * w_j
    w2 = zeros(nbmax)

    # sum of w_i  * w_j
    iw = zeros(nbmax)

    covarweight = zeros(nbmax)

    Random.seed!(iter,iseed)

    for (i, j) in iter
        # compute the distance
        for l = 1:ndims
            x0[l] = x[l][i]
            x1[l] = x[l][j]
        end
        dist = distfun(x0, x1)

        if dist > maxdist
            error("dist $(dist) is larger than maxdist $(maxdist)")
        end
        nb = floor(Int, dist / ddist + 1)
        covar[nb] =
            covar[nb] + (d[i] - datamean) * (d[j] - datamean) * weight[i] * weight[j]
        w2[nb] = w2[nb] + ((d[i] - datamean) * (d[j] - datamean))^2 * weight[i] * weight[j]
        iw[nb] = iw[nb] + weight[i] * weight[j]
    end

    covarweightmean = 0.0
    for nn = 1:nbmax
        covarweight[nn] = 0.0
        # dirty fix JMB 05/11
        # https://github.com/gher-ulg/DIVA/commit/f193cd2f5a9c350634686c730e3aa8dc606c9f59#diff-78a6698fc2fa991d95a271faa2c25d19
        if (iw[nn] > 1)
            covar[nn] = covar[nn] / iw[nn]
            w2[nn] = w2[nn] / iw[nn] - covar[nn]^2
            if (w2[nn] > 1E-8 * covar[nn]^2)
                covarweight[nn] = 1 / w2[nn]^2 * (iw[nn] - 1)
            end
            # Uniform weight
            covarweightmean = covarweightmean + covarweight[nn]
        end
    end

    for nn = 1:nbmax
        #     @show "??",nn,covarweight[nn],covarweightmean,nbmax
        covarweight[nn] = covarweightmean / nbmax + covarweight[nn]
        if (iw[nn] < 1)
            covarweight[nn] = 0
        end
    end


    # 3 iterations of a Laplacian smoother
    for jj = 1:3
        covarm = covar[1]
        for nn = 1:nbmax
            nnp = min(nbmax, nn + 1)
            covarf = covar[nn] + 0.25 * (covar[nnp] + covarm - 2 * covar[nn])
            covarm = covar[nn]
            covar[nn] = covarf
        end
    end

    ncross = 5
    RLz = -1.0
    for nn = 1:nbmax
        # if not working force simple use of variance
        if (iw[nn] != 0) && (covar[nn] < 0) && (nn > 4)
            @debug "First zero crossing: $nn $ddist $(nn*ddist)"

            RLz = ddist * nn
            ncross = nn
            break
        end
    end

    # if no zero crossing, use minimum value of covar
    if RLz == -1.0
        ncross = findmin(covar)[2]
        RLz = ddist * ncross
        @debug "No zero crossing, use minimum value at a distance of $RLz"
    end

    # Now try to fit Bessel function using only the data from ddist to zero-crossing.
    # extrapolate to zero to get S/N ratio
    @debug "Now trying to fit Bessel covariance function"
    errmin = 1.E35
    VAR = variance
    RL = RLz

    x0 = RLz / 20
    dx = ddist
    nstart = max(floor(Int, x0 / dx) + 1, 2)
    x0 = (nstart - 1) * dx
    np = floor(Int, ncross * 0.95 - 0 * nstart)
    range = nstart:(nstart+np-1)

    distx = (0:nbmax-1) * dx

    # only the distance range to be used for the optimization
    distx_range = distx[range]
    covar_range = view(covar, range)
    covarweight_range = view(covarweight, range)

    if (np < 10)
        #@show nbmax, n, nsamp, nstart, ncross
        @warn "Too few data. Will use guesses (np = $(np), RLz = $(RLz), )"
        RL = RLz
        VAR = 0.01 * variance
        SN = VAR / (variance - VAR + 1.E-10)
        varbak = 0.99 * variance
        range = nstart:(nstart+np-1)
        range = 1:0 # empty range
    else
        for ii = 1:1000
            VARtest = variance      # 17/03/2015
            RLtest = RLz / 10 + (ii - 1) * RLz / 500.0

            err, VARtest = misfit(distx_range, covar_range, covarweight_range, RLtest)

            #     @show "RL??",RLtest,VARtest,err,errmin
            if (err < errmin)
                RL = RLtest
                VAR = VARtest
                errmin = err
            end
        end

        @debug "Best fit: $RL $VAR"
        if (VAR > 0.9999 * variance)
            VAR = variance
            SN = 10000.0
            varbak = VAR
        else
            SN = VAR / (variance - VAR + 1.E-10)
        end
        @debug "S/N: $SN"
        @debug "Relative misfit of fit: $(sqrt(errmin)/VAR)"
        rqual = 1 - sqrt(errmin) / VAR
        varbak = VAR
    end

    fitcovar = VAR * (distx ./ RL) .* besselk.(1, distx ./ RL)

    dbinfo = Dict{Symbol,Any}(
        :covar => covar,
        :fitcovar => fitcovar,
        :distx => distx,
        :rqual => rqual,
        :range => range,
        :covarweight => covarweight,
        :distx => distx,
        :sn => SN,
        :meandist => meandist,
    )
    #return RL,SN,varbak,dbinfo
    stdcovar = 1 ./ sqrt.(covarweight)

    return varbak, RL, dbinfo
end

"""
    err, var = misfit(distx, covar, covarweight, RL)

This function used to be called forfit in fitlsn.f.
`err` is the weighted mean square error
`var` is the variance for a distance equal to zero.
"""
function misfit(distx, covar, covarweight, RL)
    n = length(covar)
    err = 0.0
    errb = 0.0

    # integrate the covariance and the theoretical correlation
    # over all distances
    # Their ratio is the variance

    for i = 1:n
        eps = distx[i] / RL
        errb = errb + eps * besselk(1, eps) * covarweight[i]
        err = err + covar[i] * covarweight[i]
    end

    var = err / errb

    # compute the missfit
    err = 0.0
    ww3 = 0.0

    for i = 1:n
        eps = distx[i] / RL
        covardiff = covar[i] - var * eps * besselk(1, eps)
        err = err + (covardiff^2) * covarweight[i]
        ww3 = ww3 + covarweight[i]
    end

    err = err / ww3
    return err, var
end

"""helper function for searchz"""
_getparam(z, x::Number) = x
_getparam(z, f::Function) = f(z)


"""
    lenxy,dbinfo = DIVAnd.fithorzlen(x,value,z)

Determines the horizontal correlation length `lenxy` based on the
measurements `value` at the location `x` (tuple of 3 vectors corresponding to
longitude, latitude and depth) at the depth levels defined in `z`.

Optional arguments:
 * `smoothz` (default 100): spatial filter for the correlation scale
 * `searchz` (default 50): vertical search distance (can also be a function of the depth)
 * `maxnsamp` (default 5000): maximum number of samples
 * `limitlen` (default false): limit correlation length by mean distance between
    observations
 * `limitfun` (default no function): a function with two arguments (depth and
estimated correlation length) which returns an adjusted correlation length. For
example to force the correlation length to be between 300 km and 50 km one would
use the following: `limitfun = (z,len) -> max(min(len,300_000),10_000))`. If provided
`limitfun` is used before and after the smoothing.
 * `epsilon2` (default is a vector of the same size as `value` with all elements
    equal to 1): the relative error variance of the observations. Less reliable
    observations would have a larger corresponding value.
 * `distfun`: function computing the distance between the points `xi` and `xj`.
  Per default it represents the Euclidian distance.

"""
function fithorzlen(
    x,
    value::Vector{T},
    z;
    tolrel = T(1e-4),
    smoothz = T(100.0),
    smoothk = 3,
    searchz = 50.0,
    progress = (iter, var, len, fitness) -> nothing,
    distfun = (xi, xj) -> sqrt(sum(abs2, xi - xj)),
    limitfun = (z, len) -> len,
    maxnsamp = 5000,
    limitlen = false,
    epsilon2 = ones(size(value)),
    min_rqual = 0.5,
    rng = Random.GLOBAL_RNG,
) where {T}

    if any(ϵ2 -> ϵ2 < 0, epsilon2)
        error("some values in epsilon2 are negative (minimum value is $(minimum(epsilon2)))")
    end

    kmax = length(z)
    lenopt = zeros(kmax)
    var0opt = zeros(kmax)
    fitinfos = Vector{Dict{Symbol,Any}}(undef, kmax)

    nsamp = if length(value) > maxnsamp
        maxnsamp
    else
        0 # all samples
    end

    weight = 1 ./ epsilon2
    rqual = zeros(length(z))

    Threads.@threads for k = 1:length(z)
        #for k = 1:length(z)

        searchz_k = _getparam(z[k], searchz)

        sel = if length(x) == 3
            (abs.(x[3] .- z[k]) .< searchz_k)
        else
            trues(size(x[1]))
        end

        xsel = (x[1][sel], x[2][sel])
        v = value[sel] .- mean(value[sel])

        var0opt[k], lenopt[k], fitinfos[k] =
            DIVAnd.fitlen(xsel, v, weight[sel], min(length(v), nsamp);
                          distfun = distfun,
                          rng = rng,
                          )

        rqual[k] = fitinfos[k][:rqual]
        if limitlen
            lenopt[k] = max(lenopt[k], fitinfos[k][:meandist])
        end

        @info "Data points at z=$(z[k]): $(length(v)), horz. correlation length: $(lenopt[k]) (preliminary)"
    end

    # handle layers with no data
    DIVAnd_fill!(var0opt, NaN)
    DIVAnd_fill!(lenopt, NaN)

    for k = 1:length(z)
        lenopt[k] = limitfun(z[k], lenopt[k])
    end

    # filter vertically
    lenoptf = copy(lenopt)
    rqual[rqual.<min_rqual] .= 0
    lenweight = max.(rqual, 1e-9)

    if (smoothz > 0) && (kmax > 1)
        lenoptf, lenweight = DIVAnd.smoothfilter_weighted(z, lenoptf, lenweight, smoothz)
    end
    if (smoothk > 0) && (kmax > 1)
        lenoptf, lenweight = DIVAnd.smoothfilter_weighted(z, lenoptf, lenweight, smoothk)
    end

    for k = 1:length(z)
        lenoptf[k] = limitfun(z[k], lenoptf[k])
    end

    for k = 1:length(z)
        @info "Smoothed horz. correlation length at z=$(z[k]): $(lenoptf[k])"
    end


    return lenoptf, Dict(:var0 => var0opt, :len => lenopt, :fitinfos => fitinfos)
end



"""
    lenz,dbinfo = DIVAnd.fitvertlen(x,value,z,...)

See also DIVAnd.fithorzlen
"""
function fitvertlen(
    x,
    value::Vector{T},
    z;
    smoothz = T(100.0),
    smoothk = T(3.0),
    searchz = T(10.0),
    searchxy = T(1_000.0), # meters
    maxntries = 10000,
    maxnsamp = 50,
    progress = (iter, var, len, fitness) -> nothing,
    distfun = (xi, xj) -> sqrt(sum(abs2, xi - xj)),
    limitfun = (z, len) -> len,
    epsilon2 = ones(T,size(value)),
    min_rqual = 0.5,
    rng = Random.GLOBAL_RNG,
) where {T}

    if any(ϵ2 -> ϵ2 < 0, epsilon2)
        error("some values in epsilon2 are negatives (minimum value is $(minimum(epsilon2)))")
    end

    zlevel2 = zero(T)
    zindex = Vector{Int}(undef, length(value))
    nzindex = 0

    kmax = length(z)
    lenopt = zeros(kmax)
    var0opt = zeros(kmax)
    fitinfos = Vector{Dict{Symbol,Any}}(undef, kmax)

    nsamp = min(maxnsamp, length(value))
    count = (nsamp * (nsamp - 1)) ÷ 2

    weight = 1 ./ epsilon2
    rqual = zeros(length(z))

    for k = 1:length(z)
        zlevel2 = Float64(z[k])
        searchz_k = _getparam(zlevel2, searchz)
        zindex = findall(abs.(zlevel2 .- x[3]) .< searchz_k)

        if length(zindex) == 0
            @warn "No data near z = $zlevel2"
            var0opt[k] = NaN
            lenopt[k] = NaN
            fitinfos[k] = Dict{Symbol,Any}()
        else
            iter = VertRandomCoupels(
                z[k], zindex, x,
                searchxy, maxntries, count, rng)
            #state = start(iter)
            #@code_warntype next(iter,state)
            #@code_warntype fitlen((x[3],),value,ones(size(value)),nsamp,iter)
            var0opt[k], lenopt[k], fitinfos[k] = fitlen((x[3],), value, weight, nsamp, iter; rng = rng)

            rqual[k] = fitinfos[k][:rqual]
            @info "Vert. correlation length at z=$(z[k]): $(lenopt[k])"
        end
    end

    # handle layers with no data
    DIVAnd_fill!(var0opt, NaN)
    DIVAnd_fill!(lenopt, NaN)

    for k = 1:length(z)
        lenopt[k] = limitfun(z[k], lenopt[k])
    end

    # filter vertically
    lenoptf = copy(lenopt)
    rqual[rqual.<min_rqual] .= 0
    lenweight = max.(rqual, 1e-9)

    if (smoothz > 0) && (kmax > 1)
        lenoptf, lenweight = DIVAnd.smoothfilter_weighted(z, lenoptf, lenweight, smoothz)
    end
    if (smoothk > 0) && (kmax > 1)
        lenoptf, lenweight = DIVAnd.smoothfilter_weighted(z, lenoptf, lenweight, smoothk)
    end

    for k = 1:length(z)
        lenoptf[k] = limitfun(z[k], lenoptf[k])
    end

    for k = 1:length(z)
        @debug "Smoothed vert. correlation length at z=$(z[k]): $(lenoptf[k])"
    end

    return lenoptf, Dict(:var0 => var0opt, :len => lenopt, :fitinfos => fitinfos)
end
