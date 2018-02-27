using Base.Test
import divand
using PyPlot
import WorldOceanDatabase

varname = "Salinity"
filename = "WOD-Salinity.nc"

if !isfile(filename)    
    download("https://b2drop.eudat.eu/s/UsF3RyU3xB1UM2o/download",filename)
end

value,lon,lat,depth,time,ids = divand.loadobs(Float64,filename,"Salinity")

divand.checkobs((lon,lat,depth,time),value,ids)


sel = (Dates.month.(time) .== 1)
x = (lon[sel],lat[sel],depth[sel]);
v = value[sel]
z = [0.,10,100,200,300,400,500,700,1000,1500]


srand(123);
@time lenz,infoz = divand.fitvertlen(x,v,z,len0 = 300., nmean = 500, distbin = collect([0.:50:400; 500:100:600]))


srand(123);
@time lenxy,infoxy = divand.fithorzlen(x,v,z,len0 = 3., nmean = 500, distbin = collect(0.:0.1:6))

