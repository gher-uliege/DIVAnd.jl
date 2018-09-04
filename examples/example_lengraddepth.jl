using PyPlot
using DIVAnd
using Compat: @info, range

include("./prep_dirs.jl")

dx = dy = 1.
lonr = 2.5:dx:12.
latr = 42.3:dy:44.6

bathname = joinpath(dirname(@__FILE__),"..","..","DIVAnd-example-data","Global","Bathymetry","gebco_30sec_16.nc")

bx,by,h = DIVAnd.extract_bath(bathname,true,lonr,latr)

mask = h .< 0
h[h .< 0] .= 0;

x,y = DIVAnd.ndgrid(bx,by);
pm,pn = DIVAnd.DIVAnd_metric(x,y)

L = 10_000 # m

RL = DIVAnd.lengraddepth((pm,pn),h, L)

pcolor(x,y,RL); colorbar()

figname = joinpath(figdir,basename(replace(@__FILE__,r".jl$" => ".png")));
savefig(figname)
@info "Saved figure as " * figname
