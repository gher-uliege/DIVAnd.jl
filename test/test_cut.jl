#using PyPlot
using Missings
using NCDatasets
using DataStructures
using DIVAnd

lonr = -0.2:0.1:1.1
latr = -0.2:0.1:1.1

polygon_lon = [1.05, 1.05, 0]
polygon_lat = [0, 1.05, 1.05]

sz = (length(lonr),length(latr))
field = lonr .+ latr'

filename2 = "Water_body_dissolved_oxygen_concentration_monthly2.nc"

ds = NCDataset(filename2)
lonr = ds["lon"][:]
latr = ds["lat"][:]
close(ds)

polygon_lon = [-20., 20, 20, -19]
polygon_lat = [21, 21, 53., 52]


maskkeep = DIVAnd.inpolygon(polygon_lon,polygon_lat,lonr,latr)


#filename2 = "test2.nc"
filename_cut = "cut4.nc"

#cut(filename2,varname,filename_cut,maskkeep)

DIVAnd.cut(filename2,varname,filename_cut,polygon_lon,polygon_lat)

#=
clf()
subplot(2,1,1);
#pcolormesh(lonr,latr,Float64.(maskkeep)'); colorbar()
pcolormesh(lonr,latr,field');
plot(polygon_lon[[1:end; 1]],polygon_lat[[1:end; 1]],"k-")
cl = extrema(field)
ax = axis()
clim(cl)
colorbar()

subplot(2,1,2);
pcolormesh(lonr[i],latr[j],nomissing(fieldkeep,NaN)');
plot(polygon_lon[[1:end; 1]],polygon_lat[[1:end; 1]],"k-")
clim(cl)
axis(ax)
colorbar()
=#
