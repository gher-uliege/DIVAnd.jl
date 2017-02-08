# A simple example of divand in 2 dimensions
# with observations from an analytical function.

using divand
using PyPlot

srand(1234)
# observations
nobs=99
x = rand(nobs);
y = rand(nobs);
f = sin(x*6) .* cos(y*6);
f=f+randn(nobs);

# final grid
xi,yi = ndgrid(linspace(0,1,100),linspace(0,1,100));

# reference field
fref = sin(6xi) .* cos(6yi);

# all points are valid points
mask = trues(xi);

# this problem has a simple cartesian metric
# pm is the inverse of the resolution along the 1st dimension
# pn is the inverse of the resolution along the 2nd dimension

pm = ones(xi) / (xi[2,1]-xi[1,1]);
pn = ones(xi) / (yi[1,2]-yi[1,1]);

# correlation length
len = 0.1;

# obs. error variance normalized by the background error variance
epsilon2 = 2;

# fi is the interpolated field

for imeth=0:3


bestfactorl,bestfactore, cvval,cvvalues, x2Ddata,y2Ddata,cvinter,xi2D,yi2D = divand_cv(mask,(pm,pn),(xi,yi),(x,y),f,len,epsilon2,2,3,imeth);

subplot(2,2,imeth+1)
pcolor(xi2D,yi2D,cvinter)
colorbar()
xlabel("Log10 scale factor L")
ylabel("Log10 scale factor e2")
plot(x2Ddata,y2Ddata,".")
plot(log10(bestfactorl), log10(bestfactore),"o")
title("Method $imeth")

end
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
