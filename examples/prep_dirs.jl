#=
prep_dirs.jl

Create figure and netCDF directories to store the results of the examples
Can also add other setup operations needed for to run the examples.
=#

figdir = "./figures/"
outputdir = "./netCDF/"
isdir(figdir) ? info("Directory already exists") : mkdir(figdir);
isdir(outputdir) ? info("Directory already exists") : mkdir(outputdir);
