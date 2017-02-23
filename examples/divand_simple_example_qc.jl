# A simple example of divand in 2 dimensions
# with observations from an analytical function.

using divand
using PyPlot

# observations
x = rand(300);
y = rand(300);

# Put two points in specific locations

x[1]=0.25
y[1]=0.75

x[2]=0.75
y[2]=0.25


f = sin(x*2*pi) .* sin(y*2*pi);


f=f+0.5*randn(300);

# Now fake some mix up in  two points coordinates

x[2]=0.25
y[1]=0.75

x[1]=0.75
y[2]=0.25




# final grid
xi,yi = ndgrid(linspace(0,1,30),linspace(0,1,30));


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
epsilon2 = 1;

# fi is the interpolated field
fi,s = divandrun(mask,(pm,pn),(xi,yi),(x,y),f,len,epsilon2);

qcval=divand_qc(fi,s,1)

# Find suspect points

sp=find(x-> x.>10,qcval)



suspectindexes=sortperm(qcval,rev=true)


# Copyright (C) 2014, 2017 Alexander Barth <a.barth@ulg.ac.be>
#                         Jean-Marie Beckers   <JM.Beckers@ulg.ac.be>
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
