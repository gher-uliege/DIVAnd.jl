import divand
using NCDatasets
using Interpolations
using Missings
using Base.Test

fname = "Water_body_Salinity.4Danl.nc"
varname = "Salinity"

ds = Dataset(fname)
lon = ds["lon"][:].data
lat = ds["lat"][:].data
depth = ds["depth"][:].data
time = ds["time"][:].data

v = ds["Salinity"]

i = 3
j = 3
k = 2
n = 2

loni = [lon[i]]
lati = [lat[j] ]
depthi = [depth[k]]
timei = [time[n]]


x = (lon,lat,depth)
xi = (loni,lati,depthi)

vn = zeros(size(v[:,:,:,n]))
vn[:] = map((x -> ismissing(x) ? NaN : x), v[:,:,:,n]);

fi = divand.interp(x,vn,xi)

@test fi ≈ [v[i,j,k,n]]



background = backgroundfile(fname,varname)
vn2,fi = background(xi,n,[v[i,j,k,n]],divand.Anam.notransform()[1])

@test fi ≈ [0]


