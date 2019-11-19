import DIVAnd
using Test
using DelimitedFiles
using DataStructures
using Missings
using NCDatasets

#
varname = "Salinity"
filename = "WOD-Salinity.nc"


bathname = joinpath(
    dirname(@__FILE__),
    "..",
    "..",
    "DIVAnd-example-data",
    "Global",
    "Bathymetry",
    "gebco_30sec_16.nc",
)
bathisglobal = true

obsname = joinpath(
    dirname(@__FILE__),
    "..",
    "..",
    "DIVAnd-example-data",
    "Provencal",
    "WOD-Salinity.nc",
)

cdilist = joinpath(dirname(@__FILE__), "..", "data", "CDI-list-export.csv")


if !isfile(bathname)
    @info("download bathymetry $bathname")
    bathname = download("https://dox.ulg.ac.be/index.php/s/U0pqyXhcQrXjEUX/download")
end


if !isfile(obsname)
    @info("download observations $obsname")
    obsname = download("https://dox.ulg.ac.be/index.php/s/PztJfSEnc8Cr3XN/download")
end

obsvalue, obslon, obslat, obsdepth, obstime, obsids = DIVAnd.loadobs(
    Float64,
    obsname,
    "Salinity",
)


x = (obslon,obslat)

sel  = 1:1000000

#@time weight = DIVAnd.weight_RtimesOne((obslon[sel], obslat[sel]), (0.1,0.1));
# 164.468426 seconds on dragon2 and 32 CPUs


#dup = @time DIVAnd.Quadtrees.checkduplicates((obslon[sel], obslat[sel], obsdepth[sel], obstime[sel]), obsvalue[sel], [0.1, 0.1, 0.5, 1/24], 0.01);

sel1  = 1:2000000
sel2  = (1:2000000) .+ sel1[end]


dup = @time DIVAnd.Quadtrees.checkduplicates(
    (obslon[sel1], obslat[sel1], obsdepth[sel1], obstime[sel1]), obsvalue[sel1],
    (obslon[sel2], obslat[sel2], obsdepth[sel2], obstime[sel2]), obsvalue[sel2],
    [0.1, 0.1, 0.5, 1/24], 0.01)

# 2.160539 seconds (1.03 M allocations: 2.004 GiB, 53.94% gc time)
# 1.883045 seconds (1.03 M allocations: 2.004 GiB) without GC
nothing


obsname = "/CECI/home/users/a/b/abarth/Data/Kanwal/Data_and_notebook/Global_ocean_PFL_Temperature_December.nc"

obsvalue, obslon, obslat, obsdepth, obstime, obsids = DIVAnd.loadobs(
    Float64,
    obsname,
    "Temperature",
)


sel1  = 1:2000000
sel2  = (1:2000000) .+ sel1[end]


dup = @time DIVAnd.Quadtrees.checkduplicates(
    (obslon[sel1], obslat[sel1], obsdepth[sel1], obstime[sel1]), obsvalue[sel1],
    (obslon[sel2], obslat[sel2], obsdepth[sel2], obstime[sel2]), obsvalue[sel2],
    [0.1, 0.1, 0.5, 1/24], 0.01)
