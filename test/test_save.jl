import DIVAnd
using NCDatasets
using Test
using Missings

# longitude, latitude, depth and time
xyi = (0:3, 10:13, [0, 10], [1, 2])

sz = length.(xyi)

filename = tempname()
#filename = "/tmp/test.nc"
varname = "temperature"
T = Float32
fi = randn(T, sz)
fi[1, 1, 1, 1] = NaN
mask = .!isnan.(fi)

relerr = rand(T, sz)

DIVAnd.save(filename, xyi, fi, varname; type_save = T, relerr = relerr)

ds = Dataset(filename)
fi2 = ds[varname][:,:,:,:]
close(ds)

@test fi2[mask] == fi[mask]
@test .!ismissing.(fi2) == mask

rm(filename)

# issue 68

filename = tempname()
fi = zeros(100,110)
DIVAnd.save(filename, (collect(range(0,stop=1,length=100)),
                       collect(range(0,stop=1,length=110))),
            fi, "interpolated_field")
rm(filename)
