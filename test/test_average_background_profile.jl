import DIVAnd
using Test
using DelimitedFiles
using DataStructures
using Missings
using NCDatasets
using Statistics

#
varname = "Salinity"
filename = "WOD-Salinity.nc"


obsname = joinpath(
    dirname(@__FILE__),
    "..",
    "..",
    "DIVAnd-example-data",
    "Provencal",
    "WOD-Salinity.nc",
)

if !isfile(obsname)
    @info("download observations $obsname")
    obsname = download("https://dox.ulg.ac.be/index.php/s/PztJfSEnc8Cr3XN/download")
end

obsvalue, obslon, obslat, obsdepth, obstime, obsids =
    DIVAnd.loadobs(Float64, obsname, "Salinity")



dy = dx = 0.1
lonr = 3:dx:11.8
latr = 42.0:dy:44.0
depthr = [
    0.0,
    5,
    10,
    15,
    20,
    25,
    30,
    40,
    50,
    66,
    75,
    85,
    100,
    112,
    125,
    135,
    150,
    175,
    200,
    225,
    250,
    275,
    300,
    350,
    400,
    450,
    500,
    550,
    600,
    650,
    700,
    750,
    800,
    850,
    900,
    950,
    1000,
    1050,
    1100,
    1150,
    1200,
    1250,
    1300,
    1350,
    1400,
    1450,
    1500,
    1600,
    1750,
    1850,
    2000,
];


years = 1993:1993

year_window = 10

# winter: January-March    1,2,3
# spring: April-June       4,5,6
# summer: July-September   7,8,9
# autumn: October-December 10,11,12

monthlists = [[1, 2, 3], [4, 5, 6], [7, 8, 9], [10, 11, 12]];


epsilon2 = ones(size(obsvalue))

TS = DIVAnd.TimeSelectorYW(years, year_window, monthlists)



background_filename = tempname()

if isfile(background_filename)
    rm(background_filename)
end
varname = "Salinity"


DIVAnd.average_background_profile(
    background_filename,
    (lonr, latr, depthr, TS),
    (obslon, obslat, obsdepth, obstime),
    obsvalue,
    epsilon2,
    varname,
)

@test isfile(background_filename)
nothing
