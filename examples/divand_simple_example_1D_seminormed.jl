# A simple example of divand in 1 dimensions
# with observations from an analytical function.

using divand
using PyPlot

# observations with points outside
x = collect(linspace(0,1,13))
f = sin(6*pi*x) ;

# final grid

xi=collect(linspace(-0.1,1.1,120));

# reference field
fref = sin(xi*6*pi) ;

# all points are valid points
mask = trues(size(xi));

# this problem has a simple cartesian metric
# pm is the inverse of the resolution along the 1st dimension
# pn is the inverse of the resolution along the 2nd dimension

pm = ones(size(xi)) / (xi[2]-xi[1]);


# correlation length
len = 0.1;

# obs. error variance normalized by the background error variance
epsilon2 = 1;

m = Int(ceil(1+1/2))
  # alpha is the (m+1)th row of the Pascal triangle:
  # m=0         1
  # m=1       1   1
  # m=1     1   2   1
  # m=2   1   3   3   1
  # ...
alpha = [binomial(m,k) for k = 0:m];
alpha[1]=0;

# fi is the interpolated field
fi,s = divandrun(mask,(pm,),(xi,),(x,),f,len,epsilon2;alpha=alpha);

plot(xi,fi,".",x,f,"o")
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
