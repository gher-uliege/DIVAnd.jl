# A simple example of DIVAnd in 2 dimensions
# with observations from an analytical function.

using DIVAnd
using Compat: @info, range
using PyPlot

include("./prep_dirs.jl")

# observations
x = rand(75);
y = rand(75);
f = sin.(x*6) .* cos.(y*6);

# final grid
xi,yi = ndgrid(range(0,stop=1,length=401),range(0,stop=1,length=401));

# reference field
fref = sin.(xi*6) .* cos.(yi*6);

# all points are valid points
mask = trues(size(xi));

# this problem has a simple cartesian metric
# pm is the inverse of the resolution along the 1st dimension
# pn is the inverse of the resolution along the 2nd dimension

pm = ones(size(xi)) / (xi[2,1]-xi[1,1]);
pn = ones(size(xi)) / (yi[1,2]-yi[1,1]);

# correlation length
len = 0.10;

# obs. error variance normalized by the background error variance
epsilon2 = 1.;

# fi is the interpolated field
#@time fi,s,figuess,fifine,sf,erri = DIVAndjog(mask,(pm,pn),(xi,yi),(x,y),f,len,epsilon2,ones(2),ones(2));
@time fi,s = DIVAndjog(mask,(pm,pn),(xi,yi),(x,y),f,len,epsilon2,ones(2),ones(2));

# fi is the interpolated field
@time fiex,s = DIVAndrun(mask,(pm,pn),(xi,yi),(x,y),f,len,epsilon2);

subplot(1,2,1);
pcolor(xi,yi,fi);
clim(-1,1)
colorbar()
title("result of DIVAndjog")

subplot(1,2,2);
pcolor(xi,yi,fiex);
colorbar()
clim(-1,1)
title("result of DIVAndrun")

figname = joinpath(figdir,basename(replace(@__FILE__,r".jl$" => ".png")));
savefig(figname)
@info "Saved figure as " * figname


# Copyright (C) 2014, 2018 Alexander Barth <a.barth@ulg.ac.be>
#
# This program is free software; you can redistribute it and/or modify it under
# the terms of the GNU General Public License as published by the Free Software
# Foundation; either version 2 of the License, or (at your option) any later
# version.
#
# This program is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
# FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
# details.
#
# You should have received a copy of the GNU General Public License along with
# this program; if not, see <http://www.gnu.org/licenses/>.
