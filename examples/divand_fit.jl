using Base.Test
import divand
using NLopt

min_count = 1000

fname = joinpath(dirname(@__FILE__),"..","..","divand-example-data","BlackSea","Salinity.bigfile")

value,lon,lat,depth,time,ids = divand.loadbigfile(fname)
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

figure(2)


plot(distx,covar, label = "empirical covariance");
plot(distx,fitcovar, label = "fitted function")
xlabel("normalized distance")
title("general fit")

legend()
