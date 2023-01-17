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
function fit_isotropic(
    x,
    v::Vector{T},
    distbin::Vector{T},
    mincount::Int;
    alpha = DIVAnd.alpha_default(length(x)),
    len::T = 1.0,
    minlen::T = 1e-5,
    maxlen::T = 10.0,
    tolrel::T = 1e-5,
    maxeval::Int = 10000,
    maxpoints::Int = 1000000,
    nmean::Int = 100,
    stdcovar = zeros(T, length(distbin) - 1),
    distfun = (xi, xj) -> sqrt(sum(abs2, xi - xj)),
    progress = (iter, var, len, fitness) -> nothing,
    choose = fitchoose,
) where {T}

    # number of dimensions
    n = length(x)

    # compute the empirical covariance
    # @info "Making empirical covariance"

    distx, covar, corr, varx, count, stdcovar[:] = empiriccovarmean(
        x,
        v,
        distbin,
        mincount;
        maxpoints = maxpoints,
        nmean = nmean,
        choose = choose,
        distfun = distfun,
    )

    if all(count .< mincount)
        error("Not enough pairs at all distances (count = $(count), mincount = $(mincount))")
    end
    # @info "Fitting empirical covariance"

    distx2 = copy(distx)

    # Kernel for the given dimension
    mu, K, len_scale = DIVAnd.DIVAnd_kernel(n, alpha)

    var0opt = covar[1]
    L = range(minlen, stop = maxlen, length = 10000)
    J(L) = sum(((covar - var0opt * K.(distx * len_scale / L)) ./ stdcovar) .^ 2)
    Jmin, imin = findmin(J.(L))
    lenopt = L[imin]

    fitcovar = var0opt * (K.(distx * len_scale / lenopt)::Vector{Float64})

    #return distx
    return var0opt, lenopt, distx, covar, fitcovar, stdcovar
end


function fit_isotropic(x, v, distbin, mincount; kwargs...)
    T = eltype(v)
    return fit_isotropic(x, v, collect(T, distbin), mincount; kwargs...)
end


function fitquite(iter, var0, lens, fitness) end

function fitprogress(iter, var0, lens, fitness)
    if iter == 0

        @printf("| %10s |", "var0")
        for i = 1:length(lens)
            @printf(" %10s |", "length $(i)")
        end

        @printf(" %10s |", "fitness")
        println()

        print("|------------|")
        for i = 1:length(lens)
            print("------------|")
        end

        print("------------|")
        println()
    end

    if iter % 20 == 0
        @printf("| %10g |", var0)
        for i = 1:length(lens)
            @printf(" %10g |", lens[i])
        end

        printstyled(@sprintf(" %10g ", fitness), color = :light_magenta)
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
function fit(
    x,
    v,
    distbin,
    mincount;
    alpha = DIVAnd.alpha_default(length(x)),
    minlen = zeros(length(x)),
    maxlen = ones(length(x)),
    tolrel = 1e-4,
    lens0 = ones(length(x)),
    var0 = 1.0,
    minvar0 = 0.0,
    maxvar0 = 2.0,
    maxpoints = 10000,
    nmean = 10,
    distfun = (xi, xj, lens) -> sqrt(sum(abs2, (xi - xj) ./ lens)),
    progress = (iter, var, len, fitness) -> nothing,
)
    # number of dimensions
    n = length(x)

    seed = rand(UInt64)
    iter = 0::Int

    mu, K, len_scale = DIVAnd.DIVAnd_kernel(n, alpha)

    function fitcovarlen(var0, lens)
        # declare the variable as local as they have the same name as the
        # variables in outer scope
        local distx, covar, corr, varx, count, stdcovar, fitcovar

        # fix seed to get the same observations
        Random.seed!(seed)

        distx, covar, corr, varx, count, stdcovar = empiriccovarmean(
            x,
            v,
            distbin,
            mincount;
            nmean = nmean,
            maxpoints = maxpoints,
            distfun = (xi, xj) -> distfun(xi, xj, lens),
        )

        fitcovar = var0 * K.(distx * len_scale)
        return distx, covar, fitcovar, count
    end

    function fitt(p, grad::Vector)#)
        # declare the variable as local as they have the same name as the
        # variables in outer scope

        local distx, covar, fitcovar, count, fitness
        local var0 = p[1]
        local lens = p[2:end]

        distx, covar, fitcovar, count = fitcovarlen(var0, lens)
        # avoid NaNs
        covar[isnan.(covar)] = 0
        fitness = sum(abs2, fitcovar - covar)

        progress(iter, var0, lens, fitness)
        iter = (iter::Int) + 1
        #@show var0,lens,fitness
        return fitness
    end

    n = length(x)

    #opt = Opt(:LN_COBYLA, 1+n)
    opt = Opt(:GN_DIRECT, 1 + n)
    lower_bounds!(opt, [minvar0, minlen...])
    upper_bounds!(opt, [maxvar0, maxlen...])
    xtol_rel!(opt, tolrel)

    #@show fitt([0.225,1.,1.],[])
    #@show fitt([0.225,1.,1.2],[])

    min_objective!(opt, fitt)


    minf, minx, ret = optimize(opt, [var0, lens0...])

    #@show minx
    var0opt = minx[1]
    lensopt = minx[2:end]
    distx, covar, fitcovar, count = fitcovarlen(var0opt, lensopt)

    return var0opt, lensopt, distx, covar, fitcovar
end

