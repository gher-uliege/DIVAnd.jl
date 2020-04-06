using DIVAnd
using Test

currentname = joinpath(
    dirname(@__FILE__),
    "..",
    "..",
    "DIVAnd-example-data",
    "Adriatic",
    "all-monthly-currents-2018.nc",
)

if !isfile(currentname)
    @info("download currents $currentname")
    currentname = download("https://dox.ulg.ac.be/index.php/s/qJtEotmkCZVqcx8/download")
end


yearlist = [1900:2018];
monthlist = [[1, 2, 3], [4, 5, 6], [7, 8, 9], [10, 11, 12]];
TSvelocity = DIVAnd.TimeSelectorYearListMonthList(yearlist, monthlist);

filenames = [currentname]

varnames = ("vozocrtx", "vomecrty")
outvarnames = ("u", "v")
outfilename = tempname()

DIVAnd.average_files(filenames, varnames, TSvelocity, outfilename, outvarnames)

@test isfile(outfilename)

dx, dy = 0.125, 0.125
lonr = 11.5:dx:20
latr = 39:dy:46
depthr = [0.0, 10.0, 20.0];

xi = DIVAnd.ndgrid(lonr, latr, depthr)
n = 1
veltime = [DIVAnd.ctimes(TSvelocity)[n]]

ui, vi, wi = DIVAnd.velocityfile(outfilename, ("u", "v"), TSvelocity, 1)(xi, veltime)

@test ndims(ui) == 3
