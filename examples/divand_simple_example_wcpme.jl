# A simple example of divand in 2 dimensions
# with observations from an analytical function.

using divand
using PyPlot
# using ndgrid

# observations
x = rand(75);
y = rand(75);

x = rand(10);
y = rand(10);

jmsize=120

f = sin(x*6) .* cos(y*6);

# final grid
xi,yi = ndgrid(linspace(0,1,jmsize),linspace(0,1,jmsize));

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

# obs. error variance normalized by the background error variance
epsilon2 = 0.2;
# Error scale to made comparable to the one used by divandrun in case it is not normalized
#errorscale=1 
# fi is the interpolated field
cpme=  divand_cpme(mask,(pm,pn),(xi,yi),(x,y),f,len,epsilon2);
#cpme=errorscale.*(-cpme.+1);

fi,s = divandrun(mask,(pm,pn),(xi,yi),(x,y),f,len,epsilon2);
exerr=reshape(diag(s.P),jmsize,jmsize);
# plotting of results


# plotting of results
subplot(1,2,1);
pcolor(xi,yi,cpme);
colorbar()
clim(-0.5,1.5)
plot(x,y,"k.");

subplot(1,2,2);
pcolor(xi,yi,exerr);
colorbar()
clim(-0.5,1.5)
title("Exact error");

savefig("divand_simple_example-wcpme.png")


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
