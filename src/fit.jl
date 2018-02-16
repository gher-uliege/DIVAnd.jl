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
similarily for the variable `y`. The function computes also the
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


"""
    distx,covar,corr,varx,count = empiriccovar(x,v,distbin,min_count;
                              maxpoints = 1000000,
                              distfun = (xi,xj) -> sqrt(sum(abs2,xi-xj)))

Compute the covariance, correlation and variance of a cloud of data points with
the value `v` (a vector) and the location `x` (a tuple of vectors) grouped by
distance. Random pairs are choosen and grouped by their distance
(computed by `distfun`) in bins defined by `distbin`. The function try to fill
at least `min_count` of data points in each bin but always stop after
considering `maxpoints` pairs.
"""

function empiriccovar(x,v,distbin,min_count;
                      maxpoints = 1000000,
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

    distx = (distbin[1:end-1] + distbin[2:end])/2;



    # coordinates of considered points
    xi = zeros(eltype(x[1]),length(x))
    xj = zeros(eltype(x[1]),length(x))

    for l=1:maxpoints

        # random index
        i = rand(1:length(x[1]));
        j = rand(1:length(x[1]));

        for k = 1:length(x)
            xi[k] = x[k][i]
            xj[k] = x[k][j]
        end

        if isnan(v[i]) || isnan(v[j])
            # one point is masked
            continue
        end

        #dist = distance(y(i1,j1),x(i1,j1),y(i2,j2),x(i2,j2));
        distance = distfun(xi,xj)

        if distance >= distbin[end]
            # distance too large
            continue
        end

        p = findlast(distance .>= distbin)

        if count[p] >= min_count
            # already enought points
            continue
        end

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

        if all(count .>= min_count)
            break
        end
    end

    for i = 1:pmax
        meanvi,meanvj,stdvi,stdvj,covar[i],corr[i] =
            stats(sumvi[i],sumvi2[i],sumvj[i],sumvj2[i],sumvivj[i],count[i])
        varx[i] =  ((stdvi^2  + stdvj^2) * count[i]) / (2 * count[i]-1)

        #meanx,stdx = stats(sumvi[i] + sumvj[i],sumvi2[i] + sumvj2[i],2*count[i])
        #@show i, stdx, sqrt(varx[i])
    end


    return distx,covar,corr,varx,count
end

"""
    var0,len,distx,covar,fitcovar = fit_isotropic(x,v,distbin,min_count;
                               alpha = divand.alpha_default(length(x)),
                               len = 1.,
                               var0 = 1.,
                               minlen = 0.,
                               maxlen = 10.,
                               minvar0 = 0.,
                               maxvar0 = 10.,
                               tolrel = 1e-4,
                               maxpoints = 1000000,
                               distfun = (xi,xj) -> sqrt(sum(abs2,xi-xj))),
                               progress = (var,len,fitness) -> nothing
                           )

Determines the optimal correlation length `len` and variance (for a separation
distance approaching zero) `var0` of a cloud of data points with value `v` and coordiantes
`x` (tuple of vectors with the coordinates).

The function can find the solution corresponding to  a local minimum
which is not necessarily the global minimum.

See also `empiriccovar` for future information about the output parameters.

Optional input parameters:

* `alpha`: if one correlation length is forced to zero during the anaylsis
  the values of alpha sould be set using the effective dimension.
  For example, if a 2D-analysis is simulated by forcing the vertical correlation
  length to zero, then alpha should be set to `[1,2,1]`, otherwise alpha will be
  `[1,3,3,1]` (for for any proper 3D analysis).
* `len`: initial value for the correlation length
* `var0`: initial value of the variance
* `minlen`, `maxlen`: minimum and maximum value for the correlation length
* `minvar0`, `maxvar0`: minimum and maximum value for the variance
* `tolrel`: relative tolerance for the optimizer
* `maxpoints`: maximum number of data points considered
* `distfun`: function to compute the distance between point `xi` (vector) and 
   `xj`. Per default `distun` is the Eucedian distance 
  `(xi,xj) -> sqrt(sum(abs2,xi-xj)))`.
* `progress`: call-back function to show the progress of the optimization with 
  the input parameters `var`, `len` and `fitness` (all scalars).

The length-scale parameters and the variance have the corresponding units from 
the `x` and `v`. It is therefore often necessary to provide reasonable values 
for these default parameters.

If the lower bound `minlen` is too small, then you might get the following error:

```
AmosException with id 4: input argument magnitude too large, complete loss of accuracy by argument reduction.
```

In these case, increase `minlen`.

"""
function fit_isotropic(x,v,distbin,min_count;
                       alpha = divand.alpha_default(length(x)),
                       len = 1.,
                       var0 = 1.,
                       minlen = 1e-5,
                       maxlen = 10.,
                       minvar0 = 0.,
                       maxvar0 = 10.,
                       tolrel = 1e-4,
                       maxpoints = 1000000,
                       distfun = (xi,xj) -> sqrt(sum(abs2,xi-xj)),
                       progress = (var,len,fitness) -> nothing
                       )

    # number of dimensions
    n = length(x)

    # compute the empirical covariance
    distx,covar,corr,varx,count = empiriccovar(x,v,distbin,min_count;
                                               maxpoints = maxpoints,
                                               distfun = distfun)

    # Kernel for the given dimension
    mu,K,len_scale = divand.divand_kernel(n,alpha)

    covar[count .== 0] = 0

    # function to minimize
    function fitt(p, grad::Vector #= unused =#)
        local fitcovar

        #@show p
        fitcovar = p[1] * K.(distx * len_scale/p[2])
        fitness = sum(abs2,fitcovar - covar)

        progress(p[1],p[2],fitness)
        return fitness
    end

    # setup the optimser
    opt = Opt(:LN_COBYLA, 2)
    lower_bounds!(opt, [minvar0, minlen])
    upper_bounds!(opt, [maxvar0, maxlen])
    xtol_rel!(opt,tolrel)

    min_objective!(opt, fitt)

    minf,minx,ret = optimize(opt, [var0, len])
    var0 = minx[1]
    len = minx[2]

    # fitted covariance
    fitcovar = var0 *  K.(distx * len_scale/len)

    return var0,len,distx,covar,fitcovar
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
    
    @printf("| %10g |",var0)
    for i = 1:length(lens)
        @printf(" %10g |",lens[i])
    end
                
    print_with_color(:light_magenta,@sprintf(" %10g ",fitness))
    
    println("|")
    
    return nothing
end
"""
    var0opt,lensopt,distx,covar,fitcovar = fit(x,v,distbin,min_count;
             alpha = divand.alpha_default(length(x)),
             minlen = zeros(length(x)),
             maxlen = ones(length(x)),
             tolrel = 1e-4,
             lens0 = ones(length(x)),
             var0 = 1.,
             minvar0 = 0.,
             maxvar0 = 2.,
             maxpoints = 1000000,
             distfun = (xi,xj,lens) -> sqrt(sum(abs2,(xi-xj)./lens)),
             progress = (iter,var,len,fitness) -> nothing
             )

The same as the function `fit_isotropic` except that now the correlation 
length-scale `lens0`, `minlen`, `maxlen`, `lensopt` are a vectors 
(one value per dimension). The distance function `distfun` uses an additional 
parameter to compute the normalized distance.

The note of the optional parameters in `divafit` which also applies here.
"""
function fit(x,v,distbin,min_count;
             alpha = divand.alpha_default(length(x)),
             minlen = zeros(length(x)),
             maxlen = ones(length(x)),
             tolrel = 1e-4,
             lens0 = ones(length(x)),
             var0 = 1.,
             minvar0 = 0.,
             maxvar0 = 2.,
             maxpoints = 1000000,
             distfun = (xi,xj,lens) -> sqrt(sum(abs2,(xi-xj)./lens)),
             progress = (iter,var,len,fitness) -> nothing
             )
    # number of dimensions
    n = length(x)


    const seed = rand(UInt64)
    iter = 0 :: Int
    
    mu,K,len_scale = divand.divand_kernel(n,alpha)

    function fitcovarlen(var0,lens)
        # declare the variable as local as they have the same name as the
        # variables in outer scope
        local distx, covar, corr, varx, count, fitcovar

        # fix seed to get the same observations
        srand(seed)

        distx,covar,corr,varx,count = empiriccovar(
            x,v,distbin,min_count;
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
        covar[count .== 0] = 0
        fitness = sum(abs2,fitcovar - covar)

        progress(iter,var0,lens,fitness)
        iter = (iter::Int) + 1
        #@show var0,lens,fitness
        return fitness
    end

    n = length(x)

    opt = Opt(:LN_COBYLA, 1+n)
    lower_bounds!(opt, [minvar0, minlen...])
    upper_bounds!(opt, [maxvar0, maxlen...])
    xtol_rel!(opt,tolrel)

    min_objective!(opt, fitt)

    minf,minx,ret = optimize(opt, [var0, lens0...])

    var0opt = minx[1]
    lensopt = minx[2:end]
    distx,covar,fitcovar,count = fitcovarlen(var0opt,lensopt)

    return var0opt,lensopt,distx,covar,fitcovar
end
