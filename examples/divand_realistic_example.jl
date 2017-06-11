using divand
using NetCDF

fname = joinpath(dirname(@__FILE__),"..","..","divand-example-data","BlackSea","Salinity.bigfile")
bathname = joinpath(dirname(@__FILE__),"..","..","divand-example-data","Global","Bathymetry","gebco_30sec_16.nc")


obsvalue,obslon,obslat,obsdepth,obstime,obsid = loadbigfile(fname)

dx = dy = 0.1
lonr = 27:dx:42
latr = 40:dy:47
depthr = [0., 10, 20, 30, 50, 75, 100, 125, 150, 200, 250, 300, 400, 500, 600, 700, 800, 900, 1000, 1100, 1200, 1300, 1400, 1500, 1750, 2000];
#depthr = [0, 10, 20, 30, 50, 75, 100];
#depthr = 0:10.:30.;
depthr = [0.]

timer = 1:1.:12
timer = [1]


# signal-to-noise ratio
epsilon2 = 0.1

# size of the domain

sz = (length(lonr),length(latr),length(depthr),length(timer))

# horizontal correlation length in meters
lenx = 100_000 # m
leny = 100_000 # m

# vertical correlation length in meters
lenz = Array{Float64}(sz)
for n = 1:sz[4]
    for k = 1:sz[3]
        for j = 1:sz[2]
            for i = 1:sz[1]
                lenz[i,j,k,n] = 10 + depthr[k]/5
            end
        end
    end
end

# correlation time-scale in month
lent = 1. # month

@time fi = diva(("longitude","latitude","depth","time"),
          (lonr,latr,depthr,timer),
          (obslon,obslat,obsdepth,obstime),
          obsvalue, epsilon2,
          (lenx,leny,lenz,lent),
          divand.aggregation_monthly;
          bathname = bathname
          )

#divand_save(replace(@__FILE__,r".jl$",".nc"),mask,"salinity",fi2)
nothing
