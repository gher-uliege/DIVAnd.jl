
if VERSION >= v"0.7.0-beta.0"
    using Test
    using Random
    using Dates
else
    using Base.Test
end
using DIVAnd
using Compat: @info, range
using PyPlot

varname = "Salinity"
filename = "WOD-Salinity.nc"

if !isfile(filename)
    download("https://dox.ulg.ac.be/index.php/s/PztJfSEnc8Cr3XN/download",filename)
end

value,lon,lat,depth,time,ids = DIVAnd.loadobs(Float64,filename,"Salinity")

DIVAnd.checkobs((lon,lat,depth,time),value,ids)


sel = (Dates.month.(time) .== 1)
x = (lon[sel],lat[sel],depth[sel]);
v = value[sel]
z = [0.,10,100,200,300,400,500,700,1000,1500]


if VERSION >= v"0.7.0-beta.0"
   Random.seed!(1234)
else
   srand(1234)
end
@time lenxy,infoxy = DIVAnd.fithorzlen(x,v,z)


if VERSION >= v"0.7.0-beta.0"
   Random.seed!(1234)
else
   srand(1234)
end
@time lenz,infoz = DIVAnd.fitvertlen(x,v,z; maxnsamp = 50)

nothing