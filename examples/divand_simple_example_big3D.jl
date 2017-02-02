# A simple example of divand in 3 dimensions
# with observations from an analytical function.

using divand
using PyPlot

# observations
x = rand(75);
y = rand(75);
z = rand(75);
f = sin(x*6) .* cos(y*6)+sin(z*6) .* cos(x*6) ;

# final grid
#
testsize=50
xi,yi,zi = ndgrid(linspace(0,1,testsize),linspace(0,1,testsize),linspace(0,1,testsize));

# reference field
fref = sin(xi*6) .* cos(yi*6)+sin(zi*6) .* cos(xi*6);

# all points are valid points
mask = trues(xi);

# this problem has a simple cartesian metric
# pm is the inverse of the resolution along the 1st dimension
# pn is the inverse of the resolution along the 2nd dimension

pm = ones(xi) / (xi[2,1,1]-xi[1,1,1]);
pn = ones(xi) / (yi[1,2,1]-yi[1,1,1]);
po = ones(xi) / (zi[1,1,2]-zi[1,1,1]);

# correlation length
len = 0.1;

# obs. error variance normalized by the background error variance
epsilon2 = 1;

# fi is the interpolated field
@ time fi,s = divandrun(mask,(pm,pn,po),(xi,yi,zi),(x,y,z),f,len,epsilon2);

# plotting of results
subplot(1,2,1);
pcolor(xi[:,:,15],yi[:,:,15],fref[:,:,15]);
colorbar()
clim(-1,1)
plot(x,y,"k.");

subplot(1,2,2);
pcolor(xi[:,:,15],yi[:,:,15],fi[:,:,15]);
colorbar()
clim(-1,1)
title("Interpolated field");

savefig("divand_simple_example_3D.png")

# Copyright (C) 2014, 2017 Alexander Barth <a.barth@ulg.ac.be>
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
