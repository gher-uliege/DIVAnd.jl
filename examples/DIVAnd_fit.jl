if VERSION >= v"0.7.0-beta.0"
    using Test
    using Random
    using Statistics
    using Dates
else
    using Base.Test
end
using DIVAnd
using Compat: @info, range
using PyPlot
using BenchmarkTools

mincount = 1000

fname = joinpath(dirname(@__FILE__),"..","..","DIVAnd-example-data","BlackSea","Salinity.bigfile")

value,lon,lat,depth,timed,ids = DIVAnd.loadbigfile(fname)

# 2D
# surface values for the month January
sel = (depth .> 10) .& Dates.month.(timed) .== 1;
x = (lon[sel],lat[sel]);
v = value[sel] .- mean(value[sel]);
distbin = collect(0:0.5:10)

@benchmark var0,len,distx,covar,fitcovar = DIVAnd.fit_isotropic(x,v,distbin,mincount)


#=

figure()
plot(distx,covar, label = "empirical covariance")
plot(distx,fitcovar, label = "fitted function")
legend()
title("isotropic fit 2D")


var0opt,lensopt,distx,covar,fitcovar = DIVAnd.fit(x,v,distbin,mincount)

@test all(lensopt .<= 2)

figure()
plot(distx,covar, label = "empirical covariance")
plot(distx,fitcovar, label = "fitted function")
xlabel("normalized distance")
title("general fit 2D")
legend()

# 3D
# general fit in 3D without transformation

if VERSION >= v"0.7.0-beta.0"
   Random.seed!(12345)
else
   srand(12345)
end

sel = Dates.month.(time) .== 1;


x = (lon[sel],lat[sel],depth[sel]);
v = value[sel] - mean(value[sel]);
distbin = 0:0.5:10

var0opt,lensopt,distx,covar,fitcovar =
    DIVAnd.fit(x,v,distbin,mincount;
        lens0 = [0.5, 0.5, 10.],
        maxlen = [10,10,1000],
        progress = DIVAnd.fitprogress)

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
    DIVAnd.fit(x,v,distbin,mincount;
        lens0 = [0.5, 0.5, 10.],
        maxlen = [10,10,1000],
        progress = DIVAnd.fitprogress)

figure()
plot(distx,covar, label = "empirical covariance")
plot(distx,fitcovar, label = "fitted function")
xlabel("normalized distance")
title("general fit in 3D with transformed coordinate")

println("You need to multiply the coefficient a and b of lenz by $(lensopt[3])")
=#

# Copyright (C) 2018 Jean-Marie Beckers <jm.beckers@ulg.ac.be>
#               2018 Alexander Barth <a.barth@ulg.ac.be>
#
# This program is free software; you can redistribute it and/or modify it under
# the terms of the GNU General Public License as published by the Free Software
# Foundation; either version 2 of the License, or (at your option) any later
# version.
#
# This program is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
# FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
# details.
#
# You should have received a copy of the GNU General Public License along with
# this program; if not, see <http://www.gnu.org/licenses/>.
