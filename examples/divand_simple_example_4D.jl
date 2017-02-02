# A simple example of divand in 4 dimensions
# with observations from an analytical function.

using divand
using PyPlot

# observations
nobs=200;
x = rand(nobs);
y = rand(nobs);
z = rand(nobs);
t = rand(nobs);
f = sin(x*6) .* cos(y*6)+sin(z*6) .* cos(x*6)*sin(t*2*pi) ;

# final grid
#
testsizexy=10
testsizez=5
testsizet=12
xi,yi,zi,ti = ndgrid(linspace(0,1,testsizexy),linspace(0,1,testsizexy),linspace(0,1,testsizez),linspace(0,1,testsizet));

# reference field
fref = sin(xi*6) .* cos(yi*6)+sin(zi*6) .* cos(xi*6)*sin(ti*2*pi);

# all points are valid points
mask = trues(xi);

# this problem has a simple cartesian metric
# pm is the inverse of the resolution along the 1st dimension
# pn is the inverse of the resolution along the 2nd dimension

pm = ones(xi) / (xi[2,1,1,1]-xi[1,1,1,1]);
pn = ones(xi) / (yi[1,2,1,1]-yi[1,1,1,1]);
po = ones(xi) / (zi[1,1,2,1]-zi[1,1,1,1]);
pq = ones(xi) / (ti[1,1,1,2]-ti[1,1,1,1]);

# correlation length
len = 0.1;

# obs. error variance normalized by the background error variance
epsilon2 = 1;

# fi is the interpolated field
@ time fi,s = divandrun(mask,(pm,pn,po,pq),(xi,yi,zi,ti),(x,y,z,t),f,len,epsilon2)#; moddim=[0 0 0 1]);

# plotting of results
subplot(1,2,1);
pcolor(xi[:,:,3,6],yi[:,:,3,6],fref[:,:,3,6]);
colorbar()
clim(-1,1)
plot(x,y,"k.");

subplot(1,2,2);
pcolor(xi[:,:,3,6],yi[:,:,3,6],fi[:,:,3,6]);
colorbar()
clim(-1,1)
title("Interpolated field");

savefig("divand_simple_example_$D.png")

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
