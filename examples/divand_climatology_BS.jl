#SBATCH --mem-per-cpu=8000
## ENV["MEMTOFIT"]=8.0

#=
Create a salinity climatology for the Black Sea
using different grid resolution, correlation length etc.
=#

using divand
using PyPlot
using NetCDF

include("./prep_dirs.jl");

function loadbigfile(fname)

    data = readlines(open(fname,"r"))
    nobs = length(data)

    lon = zeros(nobs)
    lat = zeros(nobs)
    depth = zeros(nobs)
    time = Array{DateTime}(nobs)
    value = zeros(nobs)
    id = Array{String}(nobs)


    for i in 1:nobs
        rec = split(data[i])
        lon[i] = parse(Float64,rec[1])
        lat[i] = parse(Float64,rec[2])
        value[i] = parse(Float64,rec[3])
        depth[i] = parse(Float64,rec[4])
        time[i] = DateTime(rec[10])
        id[i] = rec[11]
    end

    return value,lon,lat,depth,time,id
end

include("../src/override_ssmult.jl")

# if this script is in /some/path/divand.jl/examples, the data should be in
# /some/path/divand-example-data (for Linux, Mac) and likewise for Windows.
fname = joinpath(dirname(@__FILE__),"..","..","divand-example-data","BlackSea","Salinity.bigfile")
bathname = joinpath(dirname(@__FILE__),"..","..","divand-example-data","Global","Bathymetry","gebco_30sec_16.nc")
isglobal = true

if !isdefined(:value)
    value,lon,lat,depth,time,id = loadbigfile(fname)
end

@show size(value)

# Grid horizontal resolution
dx = dy = 1.
# dx = dy = 0.2
# dx = dy = 0.1
# dx = dy = 0.07
# dx = dy = 0.04
# dx = 15.0/(50*16)
# dy = 6.0/(15*16)


lonr = 27:dx:42
latr = 40.4:dy:46.6

lonr = 27:dx:42
latr = 40:dy:47


# depthr = [0., 10, 20, 30, 50, 75, 100, 125, 150, 200, 250, 300, 400, 500, 600, 700, 800, 900, 1000, 1100, 1200, 1300, 1400, 1500, 1750, 2000];
#depthr = [0, 10, 20, 30, 50, 75, 100];
depthr = 0:10.:30.;

timer = 1:1.:12

mask,(pm,pn,po,pp),(xi,yi,zi,ti) = divand_rectdom(lonr,latr,depthr,timer)
@show size(mask)     ()

epsilon2 = 1.

time2 = Dates.month.(time)

mxi,myi,mask2 = load_mask(bathname,isglobal,lonr,latr,depthr)

mask3 = repeat(mask2,inner = (1,1,1,length(timer)))


sz = size(mask)
@show sz

z = zeros(sz)
# correlation length in arc degree
#lenx = fill(1.0,sz)
#leny = fill(1.0,sz)
lenx = fill(.3,sz)
leny = fill(.3,sz)

# correlation length in meters
lenz = (10 + zi/5)/3
@show mean(lenz)
# correlation time-scale in month
lent = fill(1.,sz)

moddim=[0,0,0,12]

# Prepare background as mean vertical profile and time evolution. Just call divand in two dimensions forgetting x and y ...

vm=mean(value)
va=value-vm

toaverage=[true true false false]

@time fmb,ffb=divand_averaged_bg(mask3,(pm,pn,po,pp),(xi,yi,zi,ti),(lon,lat,depth,time2),va,(lenx,leny,4*lenz,4*lent),epsilon2*10,toaverage;moddim=moddim)


vaa=va-ffb

# Use a null correlation length along the vertical                                      v
fi,erri=divandgo(mask3,(pm,pn,po,pp),(xi,yi,zi,ti),(lon,lat,depth,time2),vaa,(lenx,leny,z,lent),epsilon2;moddim=moddim)

# Use a depth-varying correlation length along vertical                                    v
# fi,erri=divandgo(mask3,(pm,pn,po,pp),(xi,yi,zi,ti),(lon,lat,depth,time2),vaa,(lenx,leny,lenz,lent),epsilon2;moddim=moddim)


fi=fi+fmb+vm
# Why is this filter necessary; sharedArray not supported ??
#fi=divand_filter3(fi,NaN,2)
#erri=divand_filter3(erri,NaN,2)
outputfile1 = joinpath(outputdir, basename(replace(@__FILE__,r".jl$","hrS.nc")));
outputfile2 = joinpath(outputdir, basename(replace(@__FILE__,r".jl$","hrE.nc")));
divand_save(outputfile1,mask,"Salinity",fi)
divand_save(outputfile2,mask,"Errorfield",erri)
info("Results written in files\n"* outputfile1 * " and\n" * outputfile2)
nothing
