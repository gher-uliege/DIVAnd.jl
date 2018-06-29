import divand
using NCDatasets
using Interpolations
using Missings
using Base.Test

fname = "Water_body_Salinity.4Danl.nc"
varname = "Salinity"

ds = Dataset(fname)
lon = nomissing(ds["lon"][:])
lat = nomissing(ds["lat"][:])
depth = nomissing(ds["depth"][:])
time = nomissing(ds["time"][:])

v = ds["Salinity"]

i = 3
j = 2
k = 2
n = 2

loni = [lon[i]]
lati = [lat[j] ]
depthi = [10.]
timei = [time[n]]


x = (lon,lat,depth)
xi = (loni,lati,depthi)

vn = zeros(size(v[:,:,:,n]))
vn[:] = map((x -> ismissing(x) ? NaN : x), v[:,:,:,n]);

fi = divand.interp(x,vn,xi)


firef = [(v[i,j,1,n] + v[i,j,2,n])/2]
@test fi ≈ firef





background = divand.backgroundfile(fname,varname)
vn2,fi = background(xi,n,firef,divand.Anam.notransform()[1])

@test fi ≈ [0] atol=1e-5


