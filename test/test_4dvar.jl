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
mask,(pm,pn,po,pp),(xi,yi,zi,ti) = DIVAnd_squaredom(4,
    range(0,stop=1,length=7))

fi_ref = fun.(xi,yi,zi,ti)

# grid of observations
ϵ = eps()

x,y,z,t = ndgrid(range(ϵ, stop=1-ϵ, length=5),
                 range(ϵ, stop=1-ϵ, length=5),
                 range(ϵ, stop=1-ϵ, length=5),
                 range(ϵ, stop=1-ϵ, length=5));

x = x[:];
y = y[:];
z = z[:];
t = t[:];

# observations
f = fun.(x,y,z,t)

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
