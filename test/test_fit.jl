using Base.Test
import divand
using NLopt
using LsqFit

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

x = [1.,2.,3.,4.]
y = -[1.,2.,3.,4.]
sumx = sum(x)
sumx2 = sum(x.^2)

sumy = sum(y)
sumy2 = sum(y.^2)

sumxy = sum(x .* y)

meanx,stdx = stats(sumx,sumx2,length(x))
meanx,meany,stdx,stdy,covar,corr = stats(sumx,sumx2,sumy,sumy2,sumxy,length(x))

@testset "stats" begin
    @test meanx ≈ mean(x)
    @test stdx ≈ std(x)

    @test corr ≈ -1
end





function empiriccovar(x,v,distbin,min_count;
                      maxpoints = 1000000,
                      distfun = (xi,xj) -> sqrt(sum(abs2,xi-xj)))


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
    alpha: if one correlation length is forced to zero during the anaylsis
the values of alpha sould be set using the effective dimension.
For example, if a 2D-analysis is simulated by forcing the vertical correlation
length to zero, then alpha should be set to [1,2,1], otherwise alpha will be
[1,3,3,1] (for for any proper 3D analysis).

"""
function divafit(x,v,distbin,min_count;
                 alpha = divand.alpha_default(length(x)))

    len = 1.
    var0 = 1.
    @time distx,covar,corr,varx,count = empiriccovar(x,v,distbin,min_count)

    # number of dimensions
    n = length(x)

    mu,K,len_scale = divand.divand_kernel(n,alpha)
    model(x,p) = p[2] * K(x * len_scale/p[1])
    fit = curve_fit(model, distx, covar, [len,var0])

    len,var0 = fit.param

    fitcovar = var0 *  K.(distx * len_scale/len)

    return len,var0,distx,covar,fitcovar
end

"""

See the note of alpha in `divafit` which also applies here.
"""
function divafit2(x,v,distbin,min_count)
    # number of dimensions
    n = length(x)

    alpha = divand.alpha_default(length(x))
    minlen = zeros(length(x))
    maxlen = ones(length(x))
    tolrel = 1e-4
    len0 = ones(length(x))
    var0 = 1.;
    minvar0 = 0.
    maxvar0 = 2.

    const seed = rand(UInt64)

    mu,K,len_scale = divand.divand_kernel(n,alpha)

    function fitcovarlen(var0,lens)
        # fix seed to get the same observations
        srand(seed)

        distx,covar,corr,varx,count = empiriccovar(
            x,v,distbin,min_count;
            maxpoints = 1000000,
            distfun = (xi,xj) -> sqrt(sum(abs2,(xi-xj)./lens)))

        fitcovar = var0 * K.(distx * len_scale)
        return distx,covar,fitcovar,count
    end

    function fitt(p, grad::Vector #= unused =#)
        local var0 = p[1]
        local lens = p[2:3]

        distx,covar,fitcovar,count = fitcovarlen(var0,lens)
        # avoid NaNs
        covar[count .== 0] = 0
        fitness = sum(abs2,fitcovar - covar)


        @show var0,lens,fitness
        return fitness
    end

    n = length(x)

    opt = Opt(:LN_COBYLA, 1+n)
    lower_bounds!(opt, [minvar0, minlen...])
    upper_bounds!(opt, [maxvar0, maxlen...])
    xtol_rel!(opt,tolrel)

    min_objective!(opt, fitt)

    minf,minx,ret = optimize(opt, Float64[var0, len0...])
    println("got $minf at $minx after $count iterations (returned $ret)")


    var0opt = minx[1] :: Float64
    lensopt = minx[2:end] :: Vector{Float64}

    distx,covar,fitcovar,count = fitcovarlen(var0opt,lensopt)

#    distx,covar,fitcovar,count = fitcovarlen(var0,len0) :: Tuple{Vector{Float64},Vector{Float64},Vector{Float64},Vector{Int}}

    return var0opt,lensopt,distx,covar,fitcovar
end


fname = "/home/abarth/projects/Julia/divand-example-data/BlackSea/Salinity.bigfile"

value,lon,lat,depth,time,ids = divand.loadbigfile(fname)


min_count = 5000

sel = (depth .> 10) .& Dates.month.(time) .== 1;
x = (lon[sel],lat[sel]);
v = value[sel] - mean(value[sel]);
distbin = 0:0.5:10

@code_warntype divafit2(x,v,distbin,min_count)


len,var0,distx,covar,fitcovar = divafit(x,v,distbin,min_count)


plot(distx,covar, label = "empirical covariance");
plot(distx,fitcovar, label = "fitted function")
legend()


var0opt,lensopt,distx,covar,fitcovar = divafit2(x,v,distbin,min_count)

figure(2)


plot(distx,covar, label = "empirical covariance");
plot(distx,fitcovar, label = "fitted function")
xlabel("normalized distance")
legend()
