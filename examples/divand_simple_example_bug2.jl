# A simple example of divand in 2 dimensions
# with observations from an analytical function.

using divand
using PyPlot
srand(77235445)

include("./prep_dirs.jl")

# observations
x = -10+30*rand(150000);
y = -10+30*rand(150000);
f = y

# final grid
xi,yi = ndgrid(linspace(0,1,100),linspace(0,1,110));

# reference field


# all points are valid points
mask = trues(xi);

# this problem has a simple cartesian metric
# pm is the inverse of the resolution along the 1st dimension
# pn is the inverse of the resolution along the 2nd dimension

pm = ones(xi) / (xi[2,1]-xi[1,1]);
pn = ones(xi) / (yi[1,2]-yi[1,1]);

# correlation length
len = 0.03;

# obs. error variance normalized by the background error variance
epsilon2 = 0.001;

# fi is the interpolated field
@time fi,s = divandrun(mask,(pm,pn),(xi,yi),(x,y),f,len,epsilon2;alphabc=0);

# plotting of results

pcolor(xi,yi,fi);
colorbar()
clim(0,1)
#plot(x,y,"k.");
figname = joinpath(figdir,basename(replace(@__FILE__,r".jl$",".png")));
savefig(figname)
info("Saved figure as " * figname)


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
