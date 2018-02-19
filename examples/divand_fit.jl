using Base.Test
import divand
using PyPlot

mincount = 1000

fname = joinpath(dirname(@__FILE__),"..","..","divand-example-data","BlackSea","Salinity.bigfile")

value,lon,lat,depth,time,ids = divand.loadbigfile(fname)

# 2D

# surface values for the month January
sel = (depth .> 10) .& Dates.month.(time) .== 1;
x = (lon[sel],lat[sel]);
v = value[sel] - mean(value[sel]);
distbin = collect(0:0.5:10)

var0,len,distx,covar,fitcovar = divand.fit_isotropic(x,v,distbin,mincount)


figure()
plot(distx,covar, label = "empirical covariance")
plot(distx,fitcovar, label = "fitted function")
legend()
title("isotropic fit 2D")


var0opt,lensopt,distx,covar,fitcovar = divand.fit(x,v,distbin,mincount)

@test all(lensopt .<= 2)

figure()
plot(distx,covar, label = "empirical covariance")
plot(distx,fitcovar, label = "fitted function")
xlabel("normalized distance")
title("general fit 2D")
legend()

# 3D
# general fit in 3D without transformation

srand(12345)

sel = Dates.month.(time) .== 1;


x = (lon[sel],lat[sel],depth[sel]);
v = value[sel] - mean(value[sel]);
distbin = 0:0.5:10

var0opt,lensopt,distx,covar,fitcovar =
    divand.fit(x,v,distbin,mincount;
        lens0 = [0.5, 0.5, 10.],
        maxlen = [10,10,1000],
        progress = divand.fitprogress)

figure()
plot(distx,covar, label = "empirical covariance")
plot(distx,fitcovar, label = "fitted function")
xlabel("normalized distance")
title("general fit in 3D without transformation")

# with transformed coordinate

const a = 10.
const b = 1./15

# lenz and integ_lenz are related such that the
# derivative of integ_lenz is 1/lenz

lenz(z) = a + b*z
integ_lenz(z) = log.(a + b*z) / b 


x = (lon[sel],lat[sel],integ_lenz(depth[sel]));

var0opt,lensopt,distx,covar,fitcovar =
    divand.fit(x,v,distbin,mincount;
        lens0 = [0.5, 0.5, 10.],
        maxlen = [10,10,1000],
        progress = divand.fitprogress)

figure()
plot(distx,covar, label = "empirical covariance")
plot(distx,fitcovar, label = "fitted function")
xlabel("normalized distance")
title("general fit in 3D with transformed coordinate")

println("You need to multiply the coefficient a and b of lenz by $(lensopt[3])")
