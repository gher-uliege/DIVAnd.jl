using PyPlot
import divand

dx = dy = 1.
lonr = 2.5:dx:12.
latr = 42.3:dy:44.6

bathname = joinpath(dirname(@__FILE__),"..","..","divand-example-data","Global","Bathymetry","gebco_30sec_16.nc")

bx,by,b = divand.extract_bath(bathname,true,lonr,latr)

h = -b;
mask = h .< 0
h[h .< 0] = 0;

x,y = divand.ndgrid(bx,by);
pm,pn = divand.divand_metric(x,y)


L = 10_000 # m

RL = divand.lengraddepth((pm,pn),h, L)

pcolor(x,y,RL); colorbar()
