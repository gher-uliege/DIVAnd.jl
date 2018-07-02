import DIVAnd
using NCDatasets
if VERSION >= v"0.7.0-beta.0"
    using Test
else
    using Base.Test
end
using Missings

# longitude, latitude, depth and time
xyi = (0:3, 10:13, [0,10], [1,2])

sz = length.(xyi)

filename = tempname()
#filename = "/tmp/test.nc"
varname = "temperature"
T = Float32
fi = randn(T,sz)
fi[1,1,1,1] = NaN
mask = .!isnan.(fi)

relerr = rand(T,sz)

DIVAnd.save(filename,xyi,fi,varname; type_save = T, relerr = relerr)

ds = Dataset(filename)
fi2 = ds[varname][:]
close(ds)

@test fi2[mask] == fi[mask]
@test .!ismissing.(fi2) == mask
