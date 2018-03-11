# A simple example of divand in 2 dimensions
# with observations from an analytical function.

using divand
using PyPlot

include("./prep_dirs.jl")

# observations
nobs=100
x = 0.01+0.98*rand(nobs);
y = 0.01+0.98*rand(nobs);
f = sin.(x*6) .* cos.(y*6);
# x=[0.5,0.75]
# y=[0.5,0.75]
# f=[1,1]

# final grid
xi,yi = ndgrid(linspace(0,1,950),linspace(0,1,830));

# reference field
fref = sin.(xi*6) .* cos.(yi*6);

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
epsilon2 = 0.01;
vscale=0.001
vscale=0

# fi is the interpolated field
@time fiex,s = divandrun(mask,(pm,pn),(xi,yi),(x,y),f,(len,0.5*len),epsilon2;velocity=(vscale*yi,-vscale*xi),alphabc=0);

@time fi,s = divandgo(mask,(pm,pn),(xi,yi),(x,y),f,(len,0.5*len),epsilon2;velocity=(vscale*yi,-vscale*xi),alphabc=0);

subplot(1,3,1)
pcolor(xi,yi,fi)
clim(-1,1)
subplot(1,3,2)
pcolor(xi,yi,fiex)
colorbar()
clim(-1,1)
subplot(1,3,3)
pcolor(xi,yi,fiex-fi)
colorbar()

figname = joinpath(figdir,basename(replace(@__FILE__,r".jl$",".png")));
savefig(figname)
info("Saved figure as " * figname)

# Copyright (C) 2014, 2018 Alexander Barth         <a.barth@ulg.ac.be>
#                          Jean-Marie Beckers   <JM.Beckers@ulg.ac.be>
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
