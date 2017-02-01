# A simple example of divand in 2 dimensions
# with observations from an analytical function.

using divand
using PyPlot

# observations
x = -1.+3*rand(300);
y = rand(300);
f = sin(x*6) .* cos(y*6);
f=f+randn(300);

# final grid
xi,yi = ndgrid(linspace(0,1,30),linspace(0,1,30));

# reference field
fref = sin(xi*6) .* cos(yi*6);

# all points are valid points
mask = trues(size(xi));

# this problem has a simple cartesian metric
# pm is the inverse of the resolution along the 1st dimension
# pn is the inverse of the resolution along the 2nd dimension

pm = ones(size(xi)) / (xi[2,1]-xi[1,1]);
pn = ones(size(xi)) / (yi[1,2]-yi[1,1]);

# correlation length
len = 0.1;

# signal-to-noise ratio
lambda = 1;

# fi is the interpolated field
bestfact,a,b,finecv,fineloglam = divand_cvlambda(mask,(pm,pn),(xi,yi),(x,y),f,len,lambda)

plot(log10(b),a,".",fineloglam,finecv,"-",log10(bestfact),0,"o")

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
