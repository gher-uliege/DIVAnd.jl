import DIVAnd
using DelimitedFiles
using DataStructures
using Interpolations
using Missings
using NCDatasets
using Test
using Interpolations

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

obsvalue, obslon, obslat, obsdepth, obstime, obsids =
    DIVAnd.loadobs(Float64, obsname, "Salinity")

dx = dy = 0.5
lonr = 3:dx:11.8
latr = 42.0:dy:44.0
depthr = [0.0, 20.0, 30.0]
epsilon2 = 0.01


surfextend = true
sz = (length(lonr), length(latr), length(depthr))

lenx = fill(200_000, sz)
leny = fill(200_000, sz)
lenz = [50 + depthr[k] / 15 for i = 1:sz[1], j = 1:sz[2], k = 1:sz[3]]

yearlist = [1990:2000]

#monthlists = [[1, 2, 3], [4, 5, 6], [7, 8, 9], [10, 11, 12]];
monthlists = [[i] for i = 1:12]
#monthlists = [[1],[2]]
TS = DIVAnd.TimeSelectorYearListMonthList(yearlist, monthlists)

varname = "Salinity"

# File name
filename = tempname()

metadata = OrderedDict()

if isfile(filename)
    rm(filename) # delete the previous analysis
end

dbinfo = @test_logs (:info, r".*netCDF*") match_mode = :any DIVAnd.diva3d(
#dbinfo = DIVAnd.diva3d(
    (lonr, latr, depthr, TS),
    (obslon, obslat, obsdepth, obstime),
    obsvalue,
    (lenx, leny, lenz),
    epsilon2,
    filename,
    varname,
    bathname = bathname,
    bathisglobal = bathisglobal,
    surfextend = surfextend,
)


# make a analysis with the subset

index_subset = 2
monthlists_subset = [monthlists[index_subset]];
TS_subset = DIVAnd.TimeSelectorYearListMonthList(yearlist, monthlists_subset)


# File name
filename_subset = tempname()

if isfile(filename_subset)
    rm(filename_subset) # delete the previous analysis
end

dbinfo_subset = @test_logs (:info, r".*netCDF*") match_mode = :any DIVAnd.diva3d(
#dbinfo_subset = DIVAnd.diva3d(
    (lonr, latr, depthr, TS_subset),
    (obslon, obslat, obsdepth, obstime),
    obsvalue,
    (lenx, leny, lenz),
    epsilon2,
    filename_subset,
    varname,
    bathname = bathname,
    bathisglobal = bathisglobal,
    surfextend = surfextend,
)


# load data
S = NCDataset(filename) do ds
    ds[varname][:,:,:,:];
end

S_subset = NCDataset(filename_subset) do ds
    ds[varname][:,:,:,:];
end

# check if the first time instances are the same
@test all(S[:,:,:,index_subset] .=== S_subset)
