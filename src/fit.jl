"""
    meanx,stdx = stats(sumx,sumx2,N)

Computes the mean `meanx` and the standard deviation `stdx` from the sum
(`sumx`) and the sum of squares (`sumx2`) from `N` numbers.
"""
function stats(sumx,sumx2,N)
    # var(x) = 1/(N-1)  Σᵢ (xᵢ - mean(x))²
    # var(x) = 1/(N-1)  Σᵢ [ xᵢ² - 2 xᵢ mean(x) + mean(x)² ]
    # var(x) = 1/(N-1) ( Σᵢ xᵢ²  -  Σᵢ 2 xᵢ mean(x) + N mean(x)² )
    # var(x) = 1/(N-1) ( Σᵢ xᵢ²  -  2 N mean(x)² + N mean(x)² )
    # var(x) = 1/(N-1) ( Σᵢ xᵢ²  -  N mean(x)² )

    meanx = sumx/N;
    stdx = sumx2 - N * meanx.^2;

    # de-bias std
    stdx = stdx/(N-1)

    stdx =
        if stdx < 0
            zero(stdx)
        else
            sqrt(stdx)
        end

    return meanx,stdx
end

"""
    meanx,meany,stdx,stdy,covar,corr = stats(sumx,sumx2,sumy,sumy2,sumxy,N)

Computes the mean `meanx` and the standard deviation `stdx` from the sum
(`sumx`) and the sum of squares (`sumx2`) from `N` numbers and
similarly for the variable `y`. The function computes also the
Pearson correlation `corr` and covariance `covar` between `x` and `y`.
"""
function stats(sumx,sumx2,sumy,sumy2,sumxy,N)
    # (N-1) * covar(x,y) = Σᵢ (xᵢ - mean(x)) (yᵢ - mean(y))
    # (N-1) * covar(x,y) = Σᵢ (xᵢ yᵢ - xᵢ mean(y) - mean(x) yᵢ + mean(x) mean(y))
    # (N-1) * covar(x,y) = Σᵢ xᵢ yᵢ - N mean(x) mean(y) - N mean(x) mean(y) + N mean(x) mean(y)
    # (N-1) * covar(x,y) = Σᵢ xᵢ yᵢ - N mean(x) mean(y)

    meanx,stdx = stats(sumx,sumx2,N)
    meany,stdy = stats(sumy,sumy2,N)

    covar = sumxy  - sumx*sumy / N

    covar = covar/(N-1)

    corr = covar / (stdx*stdy)

    return meanx,meany,stdx,stdy,covar,corr
end


function fitchoose(x)
    # random index
    i = rand(1:length(x[1]));
    j = rand(1:length(x[1]));
    return (i,j)
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
function empiriccovar(x,v,distbin,mincount;
                      maxpoints = 10000,
                      choose = fitchoose,
                      distfun = (xi,xj) -> sqrt(sum(abs2,xi-xj)))

    @assert all(length.(x) .== length(v))


    pmax = length(distbin)-1;
    sumvivj = zeros(pmax);
    sumvi = zeros(pmax);
    sumvj = zeros(pmax);
    sumvi2 = zeros(pmax);
    sumvj2 = zeros(pmax);
    count = zeros(pmax);

    corr = zeros(pmax);
    covar = zeros(pmax);
    varx = zeros(pmax);
    count = zeros(pmax);

    distx = zeros(pmax)
    for i = 1:pmax
        distx[i] = (distbin[i] + distbin[i+1])/2
    end

    # coordinates of considered points
    xi = zeros(eltype(x[1]),length(x))
    xj = zeros(eltype(x[1]),length(x))

    for l=1:maxpoints
        # random index
        i,j = choose(x)

        for k = 1:length(x)
            xi[k] = x[k][i]
            xj[k] = x[k][j]
        end

        if isnan(v[i]) || isnan(v[j])
            # one point is masked
            continue
        end

        distance = distfun(xi,xj)

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
        sumvi[p]  += v[i]
        sumvj[p]  += v[j]
        sumvi2[p] += v[i]^2
        sumvj2[p] += v[j]^2
        count[p]  += 1;

        #if mod(l,1000) == 0
        #    @show count
        #end

        if all(count .>= mincount)
            break
        end
    end

    for i = 1:pmax
        meanvi,meanvj,stdvi,stdvj,covar[i],corr[i] =
            stats(sumvi[i],sumvi2[i],sumvj[i],sumvj2[i],sumvivj[i],count[i])
        varx[i] =  ((stdvi^2  + stdvj^2) * count[i]) / (2 * count[i]-1)

        #meanx,stdx = stats(sumvi[i] + sumvj[i],sumvi2[i] + sumvj2[i],2*count[i])
        #@show i, stdx, sqrt(varx[i])
        #@show i,distbin[i],distbin[i+1],covar[i],corr[i],sumvivj[i],count[i]
    end


    return distx,covar,corr,varx,count
end



# mean over dimensions 2 ignoring NaNs
function nm(covar::Array{T,2}) where T
    m = zeros(T,size(covar,1))
    s = zeros(T,size(covar,1))
    count = zeros(Int,size(covar,1))
    for j = 1:size(covar,2)
        for i = 1:size(covar,1)
            if !isnan(covar[i,j])
                m[i] = m[i] + covar[i,j]
                s[i] = s[i] + covar[i,j]^2
                count[i] += 1
            end
        end
    end

    for i = 1:size(covar,1)
        m[i],s[i] = stats(m[i],s[i],count[i])
    end

    return m,s
end

function empiriccovarmean(x,v::Vector{T},distbin,mincount;
                          maxpoints = 10000,
                          nmean = 10,
                          choose = fitchoose,
                          distfun = (xi,xj) -> sqrt(sum(abs2,xi-xj)))::NTuple{6,Vector{T}} where T


    sz = (length(distbin)-1,nmean)
    covar = zeros(sz)
    corr = zeros(sz)
    varx = zeros(sz)
    count = zeros(sz)

    # make sure that the variable distx is visible
    # outside the for loop
    # https://web.archive.org/save/https://docs.julialang.org/en/release-0.6/manual/variables-and-scoping/#Soft-Local-Scope-1
    # https://web.archive.org/web/20180217105449/https://stackoverflow.com/questions/22798305/function-variable-does-not-live-outside-a-for-loop

    # necessary for type inference of distx

    distx,covar[:,1],corr[:,1],varx[:,1],count[:,1] =
        empiriccovar(x,v,distbin,mincount;
                     maxpoints = maxpoints,
                     choose = choose,
                     distfun = distfun)


    for k = 2:nmean
        distx,covar[:,k],corr[:,k],varx[:,k],count[:,k] =
            empiriccovar(x,v,distbin,mincount;
                         maxpoints = maxpoints,
                         choose = choose,
                         distfun = distfun)
    end

    meancovar,stdcovar = nm(covar)
    meancorr,stdcorr = nm(corr)
    meanvarx,stdvarx = nm(varx)
    meancount,stdcount = nm(count)

    #return distx,nm(covar),nm(corr), nm(varx), nm(count)
    return distx::Vector{T},meancovar,meancorr,meanvarx,meancount,stdcovar
#    stdcovar,meancount
end

"""
    var0,len,distx,covar,fitcovar = fit_isotropic(x,v,distbin,mincount;
                               alpha = DIVAnd.alpha_default(length(x)),
                               minlen = 0.,
                               maxlen = 10.,
                               tolrel = 1e-4,
                               maxpoints = 10000,
                               nmean = 100,
                               distfun = (xi,xj) -> sqrt(sum(abs2,xi-xj))),
                               progress = (iter,var,len,fitness) -> nothing
                           )

Determines the optimal correlation length `len` and variance (for a separation
distance approaching zero) `var0` of a cloud of data points with value `v` and coordiantes
`x` (tuple of vectors with the coordinates).

The function can find the solution corresponding to a local minimum
which is not necessarily the global minimum.

See also `empiriccovar` for future information about the output parameters.

Optional input parameters:

* `alpha`: if one correlation length is forced to zero during the anaylsis
  the values of alpha sould be set using the effective dimension.
  For example, if a 2D-analysis is simulated by forcing the vertical correlation
  length to zero, then alpha should be set to `[1,2,1]`, otherwise alpha will be
  `[1,3,3,1]` (for any proper 3D analysis).
* `len`: initial value for the correlation length.
* `minlen`, `maxlen`: minimum and maximum values for the correlation length.
* `tolrel`: relative tolerance for the optimizer.
* `maxpoints`: maximum number of data points considered.
* `nmean`: the number of times an empirical covariance is estimated.
   The average covariance is used for the fitting.
* `distfun`: function to compute the distance between point `xi` (vector) and
   `xj`. Per default `distfun` is the Euclidian distance:
  `(xi,xj) -> sqrt(sum(abs2,xi-xj)))`.
* `progress`: call-back function to show the progress of the optimization with
  the input parameters `iter`, `var`, `len` and `fitness` (all scalars).

The length-scale parameters and the variance have the corresponding units from
the `x` and `v`. It is therefore often necessary to provide reasonable values
for these default parameters.

The algorithm used to estimate the correlation-length and variance is based on
randomly choosen points. Therefore the result can be different if the function
is invoked repeately. If `nmean` is increased, then these statistical
fluctuations should decrease (for a not too large value of `mincount`, i.e.
about 100 for most cases).

If the lower bound `minlen` is too small, then you might get the following error:

```
AmosException with id 4: input argument magnitude too large, complete loss of accuracy by argument reduction.
```

In these case, increase `minlen`.

"""
function fit_isotropic(x,v::Vector{T},distbin::Vector{T},mincount::Int;
                       alpha = DIVAnd.alpha_default(length(x)),
                       len::T = 1.,
                       minlen::T = 1e-5,
                       maxlen::T = 10.,
                       tolrel::T = 1e-5,
                       maxeval::Int = 10000,
                       maxpoints::Int = 1000000,
                       nmean::Int = 100,
                       stdcovar = zeros(T,length(distbin)-1),
                       distfun = (xi,xj) -> sqrt(sum(abs2,xi-xj)),
                       progress = (iter,var,len,fitness) -> nothing,
                       choose = fitchoose,
                       ) where T

    # number of dimensions
    n = length(x)

    # compute the empirical covariance
    # @info "Making empirical covariance"

    distx,covar,corr,varx,count,stdcovar[:] =
        empiriccovarmean(x,v,distbin,mincount;
                         maxpoints = maxpoints,
                         nmean = nmean,
                         choose = choose,
                         distfun = distfun)

    if all(count .< mincount)
        error("Not enough pairs at all distances (count = $(count), mincount = $(mincount))")
    end
    # @info "Fitting empirical covariance"

    distx2 = copy(distx)

    # Kernel for the given dimension
    mu,K,len_scale = DIVAnd.DIVAnd_kernel(n,alpha)

    var0opt = covar[1]
    L = range(minlen,stop=maxlen,length=10000)
    J(L) = sum(((covar - var0opt * K.(distx * len_scale/L)) ./ stdcovar).^2)
    Jmin,imin = findmin(J.(L))
    lenopt = L[imin]

    fitcovar = var0opt *  (K.(distx * len_scale/lenopt):: Vector{Float64})

    #return distx
    return var0opt,lenopt,distx,covar,fitcovar,stdcovar
end


function fit_isotropic(x,v,distbin,mincount; kwargs...)
    T = eltype(v)
    return fit_isotropic(x,v,collect(T,distbin),mincount; kwargs...)
end

function fitquite(iter,var0,lens,fitness)
end

function fitprogress(iter,var0,lens,fitness)
    if iter == 0

        @printf("| %10s |","var0")
        for i = 1:length(lens)
            @printf(" %10s |","length $(i)")
        end

        @printf(" %10s |","fitness")
        println()

        print("|------------|")
        for i = 1:length(lens)
            print("------------|")
        end

        print("------------|")
        println()
    end

    if iter % 20 == 0
        @printf("| %10g |",var0)
        for i = 1:length(lens)
            @printf(" %10g |",lens[i])
        end

        printstyled(@sprintf(" %10g ",fitness),color = :light_magenta)
        println("|")
    end

    return nothing
end
"""
    var0opt,lensopt,distx,covar,fitcovar = fit(x,v,distbin,mincount;
             alpha = DIVAnd.alpha_default(length(x)),
             minlen = zeros(length(x)),
             maxlen = ones(length(x)),
             tolrel = 1e-4,
             lens0 = ones(length(x)),
             var0 = 1.,
             minvar0 = 0.,
             maxvar0 = 2.,
             maxpoints = 10000,
             distfun = (xi,xj,lens) -> sqrt(sum(abs2,(xi-xj)./lens)),
             progress = (iter,var,len,fitness) -> nothing
             )

The same as the function `fit_isotropic` except that now the correlation
length-scale `lens0`, `minlen`, `maxlen`, `lensopt` are vectors
(one value per dimension). The distance function `distfun` uses an additional
parameter to compute the normalized distance.

The note of the optional parameters in `divafit` also applies here.
"""
function fit(x,v,distbin,mincount;
             alpha = DIVAnd.alpha_default(length(x)),
             minlen = zeros(length(x)),
             maxlen = ones(length(x)),
             tolrel = 1e-4,
             lens0 = ones(length(x)),
             var0 = 1.,
             minvar0 = 0.,
             maxvar0 = 2.,
             maxpoints = 10000,
             nmean = 10,
             distfun = (xi,xj,lens) -> sqrt(sum(abs2,(xi-xj)./lens)),
             progress = (iter,var,len,fitness) -> nothing
             )
    # number of dimensions
    n = length(x)

    seed = rand(UInt64)
    iter = 0 :: Int

    mu,K,len_scale = DIVAnd.DIVAnd_kernel(n,alpha)

    function fitcovarlen(var0,lens)
        # declare the variable as local as they have the same name as the
        # variables in outer scope
        local distx, covar, corr, varx, count, stdcovar, fitcovar

        # fix seed to get the same observations
        if VERSION >= v"0.7.0-beta.0"
   Random.seed!(seed)
else
   srand(seed)
end

        distx,covar,corr,varx,count,stdcovar =
            empiriccovarmean(
                x,v,distbin,mincount;
                nmean = nmean,
                maxpoints = maxpoints,
                distfun = (xi,xj) -> distfun(xi,xj,lens))

        fitcovar = var0 * K.(distx * len_scale)
        return distx,covar,fitcovar,count
    end

    function fitt(p, grad::Vector #= unused =#)
        # declare the variable as local as they have the same name as the
        # variables in outer scope

        local distx,covar,fitcovar,count,fitness
        local var0 = p[1]
        local lens = p[2:end]

        distx,covar,fitcovar,count = fitcovarlen(var0,lens)
        # avoid NaNs
        covar[isnan.(covar)] = 0
        fitness = sum(abs2,fitcovar - covar)

        progress(iter,var0,lens,fitness)
        iter = (iter::Int) + 1
        #@show var0,lens,fitness
        return fitness
    end

    n = length(x)

    #opt = Opt(:LN_COBYLA, 1+n)
    opt = Opt(:GN_DIRECT, 1+n)
    lower_bounds!(opt, [minvar0, minlen...])
    upper_bounds!(opt, [maxvar0, maxlen...])
    xtol_rel!(opt,tolrel)

    #@show fitt([0.225,1.,1.],[])
    #@show fitt([0.225,1.,1.2],[])

    min_objective!(opt, fitt)


    minf,minx,ret = optimize(opt, [var0, lens0...])

    #@show minx
    var0opt = minx[1]
    lensopt = minx[2:end]
    distx,covar,fitcovar,count = fitcovarlen(var0opt,lensopt)

    return var0opt,lensopt,distx,covar,fitcovar
end


function distfun_euclid(x0,x1)
    dist = zero(x0[1])
    for i = 1:length(x0)
        dist += (x0[i] - x1[i])^2
    end
    return dist = sqrt(dist)
end


distfun_m(x0,x1) = EarthRadius * distance(x0[2],x0[1],x1[2],x1[1]) * pi/180


mutable struct AllCoupels
    n:: Int
end

if VERSION >= v"0.7.0"
    function Base.iterate(iter::AllCoupels, state = (1,1))
        i,j = state
        if (i == iter.n-1) && (j == iter.n)
            return nothing
        end

        if j < iter.n
            nextstate = (i,j+1)
        else
            nextstate = (i+1,i+2)
        end

        return (nextstate,nextstate)
    end
else
    Base.start(iter::AllCoupels) = (1,1)

    function Base.next(iter::AllCoupels,state)
        i,j = state
        if j < iter.n
            nextstate = (i,j+1)
        else
            nextstate = (i+1,i+2)
        end

        return (nextstate,nextstate)
    end

    function Base.done(iter::AllCoupels,state)
        i,j = state
        return (i == iter.n-1) && (j == iter.n)
    end
end

mutable struct RandomCoupels
    n::Int
    count::Int
end

if VERSION >= v"0.7.0"
    function Base.iterate(iter::RandomCoupels, state = 0)
        if state == iter.count
            return nothing
        end

        # pick two random points
        j = rand(1:iter.n)
        i = j
        while (i == j)
            i = rand(1:iter.n)
        end

        return ((i,j),state+1)
    end
else
    Base.start(iter::RandomCoupels) = 0

    function Base.next(iter::RandomCoupels,state)
        # pick two random points
        j = rand(1:iter.n)
        i = j
        while (i == j)
            i = rand(1:iter.n)
        end
        return ((i,j),state+1)
    end
    Base.done(iter::RandomCoupels,state) = state == iter.count
end


mutable struct VertRandomCoupels
    zlevel::Float64 # depth in meters
    zindex::Vector{Int}
    x::NTuple{3,Vector{Float64}}
    searchxy::Float64 # in meters
    maxntries::Int
    count::Int
end


function _next(iter::VertRandomCoupels,state)
    # pick two random points
    j = -1
    jindex = -1

    for ntries = 1:iter.maxntries
        j = iter.zindex[rand(1:length(iter.zindex))]

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
                k = rand(1:length(iter.x[1]))
                if distfun_m([iter.x[1][k],iter.x[2][k]],[iter.x[1][j],iter.x[2][j]]) < iter.searchxy
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
        return ((j,jindex),state+1)
end

if VERSION >= v"0.7.0"
    Base.iterate(iter::VertRandomCoupels, state = 0) = (
        state == iter.count ? nothing : _next(iter,state))
else
    Base.start(iter::VertRandomCoupels) = 0
    Base.next(iter::VertRandomCoupels,state) = _next(iter,state)
    Base.done(iter::VertRandomCoupels,state) = state == iter.count
end



function fitlen(x::Tuple,d,nsamp; kwargs...)
    weight = ones(size(d))
    return fitlen(x,d,weight,nsamp; kwargs...)
end

"""
  varbak,RL,distx,covar,fitcovar,stdcovar,dbinfo =
     fitlen(x::Tuple,d,weight,nsamp; distfun = distfun_euclid, kwargs...)

this function used to be called lfit in fitlsn.f

"""
function fitlen(x::Tuple,d,weight,nsamp; kwargs...)
    # number of samples
    n = length(d)

    iter =
        if (nsamp == 0)
            AllCoupels(n)
        else
            @debug "will generate random couples"
            if (nsamp > n)
                @warn "Strange to ask for more samples than available from data; will proceed"
            end

            RandomCoupels(n,(nsamp*(nsamp-1)) ÷ 2)
        end

    if (n > 10000) && (nsamp != 0)
        @warn "Be patient big data set: $n"
    end

    return fitlen(x::Tuple,d,weight,nsamp,iter; kwargs...)
end


function fitlen(x::Tuple,d,weight,nsamp,iter; distfun = distfun_euclid, kwargs...)
    if length(d) == 0
        @warn "no data is provided to fitlen"
        return NaN,NaN,Dict{Symbol,Any}()
    end

    # number of dimensions
    ndims = length(x)

    # number of samples
    n = length(d)

    # per default operate on all data
    nop = n

    rqual = 0.
    maxdist = 0.
    meandist = 0.
    dist = 0.
    rjjj = 0.


    # compute mean and variance using the weights
    datamean = 0.
    datavar = 0.
    rn = 0.

    for i=1:n
        datamean=datamean+d[i]*weight[i]
        datavar=datavar+d[i]*d[i]*weight[i]
        rn=rn+weight[i]
    end

    datamean=datamean/rn
    variance=datavar/rn - datamean^2
    @debug "Number of data points: $n"
    @debug "data mean: $datamean"
    @debug "data variance: $variance"

    @debug "Now calculating distance distribution"

    x0 = zeros(ndims)
    x1 = zeros(ndims)

    if VERSION >= v"0.7.0-beta.0"
        Random.seed!(n)
    else
        srand(n)
    end

    for (i,j) in iter
        # compute the distance
        for l = 1:ndims
            x0[l] = x[l][i]
            x1[l] = x[l][j]
        end
        dist = distfun(x0,x1)

        meandist = meandist+dist
        if (dist > maxdist)
            maxdist = dist
        end
    end

    if (nsamp == 0)
        rjjj = rn*(rn-1.)*0.5
    else
        rjjj = nsamp*(nsamp-1) / 2
    end

    @debug "Number of data couples considered: $rjjj"
    meandist=meandist/rjjj

    @debug "maximum distance between points: $maxdist"

    @debug "Mean distance between points: $meandist"

    rnbins =
        if (nsamp == 0)
            min(80.,rn^2 / maxdist * meandist / 20.)
        else
            min(80.,nsamp^2 / maxdist * meandist / 20.)
        end

    @debug "Number of probable active bins: $rnbins"

    ddist = meandist/rnbins
    nbmax = floor(Int,maxdist / ddist + 1)
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

    if VERSION >= v"0.7.0-beta.0"
        Random.seed!(n)
    else
        srand(n)
    end

    for (i,j) in iter
        # compute the distance
        for l = 1:ndims
            x0[l] = x[l][i]
            x1[l] = x[l][j]
        end
        dist = distfun(x0,x1)

        if dist > maxdist
            error("dist $(dist) is larger than maxdist $(maxdist)")
        end
        nb = floor(Int,dist/ddist+1)
        covar[nb] = covar[nb] + (d[i]-datamean)*(d[j]-datamean) * weight[i]*weight[j]
        w2[nb] = w2[nb] +((d[i]-datamean)*(d[j]-datamean))^2 * weight[i]*weight[j]
        iw[nb] = iw[nb] + weight[i]*weight[j]
    end

    covarweightmean = 0.
    for nn=1:nbmax
        covarweight[nn] = 0.
        # dirty fix JMB 05/11
        # https://github.com/gher-ulg/DIVA/commit/f193cd2f5a9c350634686c730e3aa8dc606c9f59#diff-78a6698fc2fa991d95a271faa2c25d19
        if (iw[nn] > 1)
            covar[nn]=covar[nn]/iw[nn]
            w2[nn]=w2[nn]/iw[nn]-covar[nn]^2
            if (w2[nn] > 1E-8*covar[nn]^2)
                covarweight[nn]=1/w2[nn]^2*(iw[nn]-1)
            end
            # Uniform weight
            covarweightmean=covarweightmean+covarweight[nn]
        end
    end

    for nn=1:nbmax
        #     @show "??",nn,covarweight[nn],covarweightmean,nbmax
        covarweight[nn]=covarweightmean/nbmax+covarweight[nn]
        if (iw[nn] < 1)
            covarweight[nn]=0
        end
    end


    # 3 iterations of a Laplacian smoother
    for jj=1:3
        covarm=covar[1]
        for nn=1:nbmax
            nnp=min(nbmax,nn+1)
            covarf=covar[nn]+0.25*(covar[nnp]+covarm-2*covar[nn])
            covarm=covar[nn]
            covar[nn]=covarf
        end
    end

    ncross = 5
    RLz = -1.
    for nn=1:nbmax
        # if not working force simple use of variance
        if (iw[nn] != 0) && (covar[nn] < 0) && (nn > 4)
            @debug "First zero crossing: $nn $ddist $(nn*ddist)"

            RLz = ddist*nn
            ncross = nn
            break
        end
    end

    # if no zero crossing, use minimum value of covar
    if RLz == -1.
        ncross = findmin(covar)[2]
        RLz = ddist*ncross
        @debug "No zero crossing, use minimum value at a distance of $RLz"
    end

    # Now try to fit Bessel function using only the data from ddist to zero-crossing.
    # extrapolate to zero to get S/N ratio
    @debug "Now trying to fit Bessel covariance function"
    errmin = 1.E35
    VAR = variance
    RL = RLz

    x0 = RLz/20
    dx = ddist
    nstart = max(floor(Int,x0 / dx)+1,2)
    x0 = (nstart-1)*dx
    np = floor(Int,ncross*0.95-0*nstart)
    range = nstart:(nstart+np-1)

    distx = (0:nbmax-1) * dx

   # only the distance range to be used for the optimization
   distx_range = distx[range]
   covar_range = view(covar,range)
   covarweight_range = view(covarweight,range)

    if (np < 10)
        #@show nbmax, n, nsamp, nstart, ncross
        @warn "Too few data. Will use guesses (np = $(np), RLz = $(RLz), )"
        RL=RLz
        VAR=0.01*variance
        SN=VAR/(variance-VAR+1.E-10)
        varbak = 0.99*variance
        range = nstart:(nstart+np-1)
        range = 1:0 # empty range
    else
        for ii=1:1000
            VARtest=variance      # 17/03/2015
            RLtest=RLz/10+(ii-1)*RLz/500.

            err,VARtest = misfit(distx_range,covar_range,covarweight_range,RLtest)

            #     @show "RL??",RLtest,VARtest,err,errmin
            if (err < errmin)
                RL=RLtest
                VAR=VARtest
                errmin=err
            end
        end

        @debug "Best fit: $RL $VAR"
        if (VAR > 0.9999*variance)
            VAR = variance
            SN = 10000.
            varbak = VAR
        else
            SN=VAR/(variance-VAR+1.E-10)
        end
        @debug "S/N: $SN"
        @debug "Relative misfit of fit: $(sqrt(errmin)/VAR)"
        rqual=1-sqrt(errmin)/VAR
        varbak=VAR
    end

    fitcovar = VAR * (distx./RL) .* besselk.(1,distx./RL)

    dbinfo = Dict{Symbol,Any}(
        :covar => covar,
        :fitcovar => fitcovar,
        :distx => distx,
        :rqual => rqual,
        :range => range,
        :covarweight => covarweight,
        :distx => distx,
        :sn => SN,
        :meandist => meandist
    )
    #return RL,SN,varbak,dbinfo
    stdcovar = 1 ./ sqrt.(covarweight)

    return varbak,RL,dbinfo
end

# this function used to be called forfit in fitlsn.f

function misfit(distx,covar,covarweight,RL)
    n = length(covar)
    err = 0.
    errb = 0.

    # integrate the covariance and the theoretical correlation
    # over all distances
    # Their ratio is the variance

    for i=1:n
        eps = distx[i]/RL
        errb = errb + eps*besselk(1,eps) * covarweight[i]
        err = err + covar[i]*covarweight[i]
    end

    var=err/errb

    # compute the missfit
    err = 0.
    ww3=0.

    for i=1:n
        eps = distx[i]/RL
        covardiff = covar[i] - var*eps*besselk(1,eps)
        err = err + (covardiff^2) * covarweight[i]
        ww3 = ww3 + covarweight[i]
      end

    err=err/ww3
    return err,var
end



"""
    lenz,dbinfo = DIVAnd.fithorzlen(x,value,z)

Determines the horizontal correlation length `lenz` based on the
measurments `value` at the location `x` (tuple of 3 vectors corresponding to
longitude, latitude and depth).

Optional arguments:
 * `smoothz` (default 100): spatial filter for the correlation scale
 * `searchz` (default 50): vertical search distance
 * `maxnsamp` (default 5000): maximum number of samples
 * `limitlen` (default false): limit correlation length by mean distance between
    observations
 * `limitfun` (default no function): a function with with the two arguments (depth and
estimated correlation length) which returns an adjusted correlation length. For
example to force the correlation length to be between 300 km and 50 km one would
use the following: `limitfun = (z,len) -> max(min(len,300),10))`. If provided
`limitfun` is used before and after the smoothing.
 * `epsilon2` (default is a vector of the same size as `value` with all elements
    equal to 1): the relative error variance of the observations. Less reliable
    observation would have a larger corresponding value.


"""
function fithorzlen(x,value::Vector{T},z;
                    tolrel::T = 1e-4,
                    smoothz::T = 100.,
                    smoothk::T = 3.,
                    searchz::T = 50.,
                    progress = (iter,var,len,fitness) -> nothing,
                    distfun = (xi,xj) -> sqrt(sum(abs2,xi-xj)),
                    limitfun = (z,len) -> len,
                    maxnsamp = 5000,
                    limitlen = false,
                    epsilon2 = ones(size(value)),
                    ) where T

    if any(ϵ2 -> ϵ2 < 0,epsilon2)
        error("some values in epsilon2 are negatives (minimum value is $(minimum(epsilon2)))")
    end

    kmax = length(z)
    lenopt = zeros(kmax)
    var0opt = zeros(kmax)
    fitinfos = Vector{Dict{Symbol,Any}}(undef,kmax)

    nsamp =
        if length(value) > maxnsamp
            maxnsamp
        else
            0 # all samples
        end

    weight = 1 ./ epsilon2

    for k = 1:length(z)

        sel =
            if length(x) == 3
                (abs.(x[3] .- z[k]) .< searchz)
            else
                trues(size(x[1]))
            end

        xsel = (x[1][sel],x[2][sel]);
        v = value[sel] .- mean(value[sel]);

        var0opt[k],lenopt[k],fitinfos[k] = DIVAnd.fitlen(
            xsel,v,weight[sel],nsamp;
            distfun = distfun
        )

        if limitlen
            lenopt[k] = max(lenopt[k], fitinfos[k][:meandist])
        end

        @info "Data points at z=$(z[k]): $(length(v)), horz. correlation length: $(lenopt[k])"
    end

    # handle layers with no data
    DIVAnd_fill!(var0opt,NaN)
    DIVAnd_fill!(lenopt,NaN)

    for k = 1:length(z)
        lenopt[k] = limitfun(z[k],lenopt[k])
    end

    # filter vertically
    lenoptf = copy(lenopt)
    if (smoothz > 0) && (kmax > 1)
        DIVAnd.smoothfilter!(z,lenoptf,smoothz)
    end
    if (smoothk > 0) && (kmax > 1)
        DIVAnd.smoothfilter!(1:length(z),lenoptf,smoothk)
    end

    for k = 1:length(z)
        lenoptf[k] = limitfun(z[k],lenoptf[k])
    end

    for k = 1:length(z)
        @debug "Smoothed horz. correlation length at z=$(z[k]): $(lenoptf[k])"
    end


    @show z,lenoptf

    return lenoptf,Dict(
        :var0 => var0opt,
        :len => lenopt,
        :fitinfos => fitinfos,
   )
end



"""
    lenz,dbinfo = DIVAnd.fitvertlen(x,value,z,...)

See also DIVAnd.fithorzlen
"""
function fitvertlen(x,value::Vector{T},z;
                     smoothz::T = 100.,
                     smoothk::T = 3.,
                     searchz::T = 10.,
                     searchxy::T = 1_000., # meters
                     maxntries::Int = 10000,
                     maxnsamp = 50,
                     progress = (iter,var,len,fitness) -> nothing,
                     distfun = (xi,xj) -> sqrt(sum(abs2,xi-xj)),
                     limitfun = (z,len) -> len,
                     epsilon2 = ones(size(value)),
                    ) where T

    if any(ϵ2 -> ϵ2 < 0,epsilon2)
        error("some values in epsilon2 are negatives (minimum value is $(minimum(epsilon2)))")
    end

    zlevel2 = zero(T)
    zindex = Vector{Int}(undef,length(value))
    nzindex = 0

    kmax = length(z)
    lenopt = zeros(kmax)
    var0opt = zeros(kmax)
    fitinfos = Vector{Dict{Symbol,Any}}(undef,kmax)

    nsamp = min(maxnsamp,length(value))
    count = (nsamp*(nsamp-1)) ÷ 2

    weight = 1 ./ epsilon2

    for k = 1:length(z)
        zlevel2 = Float64(z[k])
        zindex = findall(abs.(zlevel2 .- x[3]) .< searchz)

        if length(zindex) == 0
            @warn "No data near z = $zlevel2"
            var0opt[k] = NaN
            lenopt[k] = NaN
            fitinfos[k] = Dict{Symbol,Any}()
        else
            iter = VertRandomCoupels(z[k],zindex,x,searchxy,maxntries,count)
            #state = start(iter)
            #@code_warntype next(iter,state)
            #@code_warntype fitlen((x[3],),value,ones(size(value)),nsamp,iter)
            var0opt[k],lenopt[k],fitinfos[k] = fitlen((x[3],),value,weight,nsamp,iter)

            @info "Vert. correlation length at z=$(z[k]): $(lenopt[k])"
        end
    end

    # handle layers with no data
    DIVAnd_fill!(var0opt,NaN)
    DIVAnd_fill!(lenopt,NaN)

    for k = 1:length(z)
        lenopt[k] = limitfun(z[k],lenopt[k])
    end

    # filter vertically
    lenoptf = copy(lenopt)
    if (smoothz > 0) && (kmax > 1)
        DIVAnd.smoothfilter!(z,lenoptf,smoothz)
    end
    if (smoothk > 0) && (kmax > 1)
        DIVAnd.smoothfilter!(1:length(z),lenoptf,smoothk)
    end

    for k = 1:length(z)
        lenoptf[k] = limitfun(z[k],lenoptf[k])
    end

    for k = 1:length(z)
        @debug "Smoothed vert. correlation length at z=$(z[k]): $(lenoptf[k])"
    end

    return lenoptf,Dict(
        :var0 => var0opt,
        :len => lenopt,
        :fitinfos => fitinfos,
   )

end
