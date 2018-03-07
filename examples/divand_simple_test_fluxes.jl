# A simple example of divand in 2 dimensions
# with observations from an analytical function.

using divand
using PyPlot

# function to interpolate
fun(x,y) = sin.(6x) * cos.(6y)

# observations

x = rand(1);
y = rand(1);
f = fun.(x,y)

# final grid
xi,yi = ndgrid(linspace(0,100,100),linspace(0,110,110));

# reference field
fref = fun.(xi,yi)

# all points are valid points
mask = trues(xi);

mask[1,:]=false;
mask[end,:]=false;

mask[30:80,30:80]=false;

# this problem has a simple cartesian metric
# pm is the inverse of the resolution along the 1st dimension
# pn is the inverse of the resolution along the 2nd dimension

pm = ones(xi) / (xi[2,1]-xi[1,1]);
pn = ones(xi) / (yi[1,2]-yi[1,1]);

# correlation length
len = 10;

# obs. error variance normalized by the background error variance
epsilon2 = 10000000.;

h=xi.*(100.-xi)+20

# fi is the interpolated field

fluxes=sin.(yi[1,:]./10.)+0.1*rand(size(h)[2])

@time fi,s = divandrun(mask,(pm,pn),(xi,yi),(x,y),f,len,epsilon2;topographyforfluxes=h,fluxes=fluxes,epsfluxes=1.,alphabc=1,alpha=[1, 0, 1]);

# plotting of results

fluxesafter=zeros(size(h)[2])

for j=1:size(h)[2]
 for i=2:size(h)[1]-2
	if mask[i,j]&& mask[i+1,j]
 		fluxesafter[j]=fluxesafter[j]+h[i,j]*(fi[i+1,j]-fi[i,j])
	end
 end
end
 
 

pcolor(xi,yi,fi);
colorbar()

title("Interpolated field");

savefig("divand_simple_example_fluxes.png")

var(fluxes-fluxesafter)

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
