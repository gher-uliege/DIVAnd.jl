# Testing DIVAnd in 2 dimensions with advection.
using DIVAnd
if VERSION >= v"0.7.0-beta.0"
    using Test
else
    using Base.Test
end

# grid of background field
mask,(pm,pn),(xi,yi) = DIVAnd_squaredom(
    2,range(-1, stop = 1, length = 30))

# island at these location
mi0 = 12
mi1 = 13

mj0 = 12
mj1 = 13

mask[mi0:mi1,mj0:mj1] .= false

x = [0.]
y = [0.]
f = [1.]

epsilon2 = 1/200
len = 0.3

fi,s = DIVAndrun(mask,(pm,pn),(xi,yi),(x,y),f,len,epsilon2);

eps2 = 1e-7

c = DIVAnd_constr_constcoast(mask,eps2)

@test DIVAnd.gradcoast(mask,fi[mask]) ≈ DIVAnd.sparse_gradcoast(mask) * fi[mask]


fi2,s = DIVAndrun(mask,(pm,pn),(xi,yi),(x,y),f,len,epsilon2;
                  constraints = [c]);


# extract boundary values around the island and check
# if the standard deviation is small

fi_coast = fi2[mi0-1:mi1+1,mj0-1:mj1+1][mask[mi0-1:mi1+1,mj0-1:mj1+1]]
@test fi_coast ≈ fill(fi_coast[1],size(fi_coast)) atol=1e-5


# more complex example

if VERSION >= v"0.7.0-beta.0"
   Random.seed!(1234)
else
   srand(1234)
end
mask0,(pm,pn),(xi,yi) = DIVAnd_squaredom(
    2,range(0, stop = 1, length = 100))
mask = DIVAnd.random(mask0,(pm,pn),0.1,1)[:,:,1] .> 0.5
x = rand(Float64,100)
y = rand(Float64,size(x))
f = sin.(2*π*x) .* sin.(2*π*y)

len = 0.1

fi3,s = DIVAndrun(mask,(pm,pn),(xi,yi),(x,y),f,len,epsilon2;
                  constraints = [DIVAnd_constr_constcoast(mask,eps2)]);


nothing



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
