if VERSION >= v"0.7.0-beta.0"
    using Test
else
    using Base.Test
end

using DIVAnd

bathname = joinpath(dirname(@__FILE__),"..","..","DIVAnd-example-data","Global","Bathymetry","gebco_30sec_16.nc")

x0 = 27.0
x1 = 41.8
dx = 0.12
y0 = 40.3
y1 = 46.8
dy = 0.1
level  = 0
isglobal = true;


xi,yi,mi = load_mask(bathname,isglobal,x0,x1,dx,y0,y1,dy,level)


levels = [0,100,1000]
xi,yi,mis = load_mask(bathname,isglobal,x0,x1,dx,y0,y1,dy,levels)

@test mi == mis[:,:,1]
