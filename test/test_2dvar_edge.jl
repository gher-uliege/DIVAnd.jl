# Testing DIVAnd in 2 dimensions with independent verification.

if VERSION >= v"0.7.0-beta.0"
    using Test
else
    using Base.Test
end
using DIVAnd

scalefactor = 3

# grid of background field
#mask,(pm,pn),(xi,yi) = DIVAnd_squaredom(2,scalefactor*range(0,stop=1,length=4),scalefactor*range(0,stop=1,length=))
mask,(pm,pn),(xi,yi) = DIVAnd_rectdom(scalefactor*range(0,stop=1,length=4),scalefactor*range(0,stop=1,length=3))
mask[:,end] .= false
mask[2:2,end] .= true


epsilon = 1e-10;

# grid of observations

x,y = ndgrid(scalefactor * range(0.1,stop = 0.9, length = 2),
             scalefactor * [0.3])

v = copy(x)/maximum(x)
#v[:,2] = reverse(v[:,2])

x = x[:]
y = y[:]
v = v[:]

lenx = scalefactor * 100.;
leny = scalefactor * 100.;

epsilon2 = 0.0001;

#,err,s
va,s = DIVAndrun(mask,(pm,pn),(xi,yi),(x,y),v,(lenx,leny),epsilon2,primal=true)

@test any(abs.(va[isfinite.(va)]) .> 1)


va,s = DIVAndrun(mask,(pm,pn),(xi,yi),(x,y),v,(lenx,leny),epsilon2,primal=true,alphabc=0,coeff_derivative2 = [1.,1.]);
@test all(abs.(va[isfinite.(va)]) .< 1)


# Copyright (C) 2014, 2016 Alexander Barth <a.barth@ulg.ac.be>
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
