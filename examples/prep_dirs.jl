#=
prep_dirs.jl

Create figure and netCDF directories to store the results of the examples
Can also add other setup operations needed for to run the examples.
=#

figdir = "./figures/"
outputdir = "./netCDF/"

if isdir(figdir)
    @info "Directory already exists"
else
    mkdir(figdir)
end

if isdir(outputdir)
    @info "Directory already exists"
else
    mkdir(outputdir)
end
