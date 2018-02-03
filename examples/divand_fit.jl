using Base.Test
import divand
using NLopt

include("../src/fit.jl")
min_count = 1000



fname = joinpath(dirname(@__FILE__),"..","..","divand-example-data","BlackSea","Salinity.bigfile")

value,lon,lat,depth,time,ids = divand.loadbigfile(fname)

#=
# 2D

sel = (depth .> 10) .& Dates.month.(time) .== 1;
x = (lon[sel],lat[sel]);
v = value[sel] - mean(value[sel]);
distbin = 0:0.5:10


#@code_warntype fit(x,v,distbin,min_count)


len,var0,distx,covar,fitcovar = fit_isotropic(x,v,distbin,min_count)


figure(1)
plot(distx,covar, label = "empirical covariance");
plot(distx,fitcovar, label = "fitted function")
legend()
title("isotropic fit")


var0opt,lensopt,distx,covar,fitcovar = fit(x,v,distbin,min_count)

@test all(lensopt .<= 2)

figure()
plot(distx,covar, label = "empirical covariance");
plot(distx,fitcovar, label = "fitted function")
xlabel("normalized distance")
title("general fit")

legend()
=#

# 3D

srand(12345)

sel = Dates.month.(time) .== 1;
const a = 10.
const b = 1./15

#ll = 12.3
#const a = 10. * 12.3
#const b = 1./15  * 12.3

# lenz and integ_lenz are related such that the
# derivative of integ_lenz is 1/lenz

lenz(z) = a + b*z
integ_lenz(z) = log.(a + b*z) / b 


x = (lon[sel],lat[sel],depth[sel]);
v = value[sel] - mean(value[sel]);
distbin = 0:0.5:10

#=
var0opt,lensopt,distx,covar,fitcovar =
    fit(x,v,distbin,min_count;
        lens0 = [0.5, 0.5, 10.],
        maxlen = [10,10,1000],
        progress = fitprogress)

figure()
plot(distx,covar, label = "empirical covariance");
plot(distx,fitcovar, label = "fitted function")
xlabel("normalized distance")
title("general fit in 3D")
=#

# with transformed coordinate

x = (lon[sel],lat[sel],integ_lenz(depth[sel]));

#        lens0 = [0.5, 0.5, 10. ./ ll],
#        maxlen = [10,10,1000 ./ ll],

var0opt,lensopt,distx,covar,fitcovar =
    fit(x,v,distbin,min_count;
        lens0 = [0.5, 0.5, 10.],
        maxlen = [10,10,1000],
        progress = fitprogress)

figure()
plot(distx,covar, label = "empirical covariance");
plot(distx,fitcovar, label = "fitted function")
xlabel("normalized distance")
title("general fit in 3D with transformed coordinate")

println("You need to multiply the coefficient a and b of lenz by $(lensopt[3])")
