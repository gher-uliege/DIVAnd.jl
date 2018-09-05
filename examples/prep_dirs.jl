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



# work-around for PyPlot
# https://github.com/JuliaPy/PyPlot.jl/issues/362
if VERSION >= v"0.7"
    import LinearAlgebra
    import PyPlot: pcolor

    pcolor(x::AbstractRange, y::AbstractRange, z::LinearAlgebra.Adjoint{Float64,Array{Float64,2}}; kwargs...) = pcolor(x,y,copy(z); kwargs...)
    pcolor(x,y,z::LinearAlgebra.Adjoint{Float64,Array{Float64,2}}; kwargs...) = pcolor(x,y,copy(z); kwargs...)
    pcolor(z::LinearAlgebra.Adjoint{Float64,Array{Float64,2}}; kwargs...) = pcolor(copy(z); kwargs...)
end