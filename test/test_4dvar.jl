# Testing DIVAnd in 4 dimensions.

if VERSION >= v"0.7.0-beta.0"
    using Test
else
    using Base.Test
end

# function to interpolate
k = π
fun(x,y,z,t) = sin(k*x) * sin(k*y) * sin(k*z) * sin(k*t)

# grid of background field
xi,yi,zi,ti = ndgrid(linspace(0,1.,7),linspace(0,1.,7),linspace(0,1.,7),linspace(0,1.,7));
fi_ref = fun.(xi,yi,zi,ti)

# grid of observations
ϵ = eps()
x,y,z,t = ndgrid(linspace(ϵ,1-ϵ,5),linspace(ϵ,1-ϵ,5),linspace(ϵ,1-ϵ,5),linspace(ϵ,1-ϵ,5));
x = x[:];
y = y[:];
z = z[:];
t = t[:];

# observations
f = fun.(x,y,z,t)

# all points are valid points
mask = trues(xi);

# this problem has a simple cartesian metric
# pm (pn,po,pp) is the inverse of the resolution along the 1st (2nd, 3rd, 4th) dimension
pm = ones(xi) / (xi[2,1,1,1]-xi[1,1,1,1]);
pn = ones(xi) / (yi[1,2,1,1]-yi[1,1,1,1]);
po = ones(xi) / (zi[1,1,2,1]-zi[1,1,1,1]);
pp = ones(xi) / (ti[1,1,1,2]-zi[1,1,1,1]);

# correlation length
len = 0.12;

# obs. error variance normalized by the background error variance
epsilon2 = 0.1;

# fi is the interpolated field
fi,s = DIVAndrun(mask,(pm,pn,po,pp),(xi,yi,zi,ti),(x,y,z,t),f,len,epsilon2);

# compute RMS to background field
rms = sqrt(mean((fi_ref[:] - fi[:]).^2));

# rms should be 0.0111807

@test rms < 0.012


# Copyright (C) 2014,2017 Alexander Barth <a.barth@ulg.ac.be>
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
