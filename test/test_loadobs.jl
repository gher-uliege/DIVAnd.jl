if VERSION >= v"0.7.0-beta.0"
    using Test
    using Dates
    using Random
else
    using Base.Test
end
import DIVAnd


T = Float32

fname = tempname() * ".nc"
nobs = 10
lon = randn(T,nobs)
lat = randn(T,nobs)
depth = randn(T,nobs)
time = rand(DateTime(2000,1,1):Dates.Day(1):DateTime(2010,1,1),nobs)
value = randn(T,nobs)
ids = String[randstring(10) for i in 1:nobs]


DIVAnd.saveobs(fname,"Salinity",value,(lon,lat,depth,time),ids;type_save = T)

value2,lon2,lat2,depth2,time2,ids2 = DIVAnd.loadobs(T,fname,"Salinity")

@test value == value2
@test lon == lon2
@test lat == lat2
@test depth == depth2
@test ids == ids2

rm(fname)

