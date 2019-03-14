# Testing DIVAnd in 2 dimensions with independent verification.

if VERSION >= v"0.7.0-beta.0"
    using Test
else
    using Base.Test
end
using DIVAnd


# grid of background field
mask,(pm,pn),(xi,yi) = DIVAnd_squaredom(2,0.:1:2)
mask .= false
mask[2,2] = true

epsilon = 1e-10;

# grid of observations
nobs = 10000
x = fill(1.99,(nobs,))
y = fill(1.,(nobs,))
v = fill(1.,(nobs,))

lenx = .15;
leny = .15;

epsilon2 = 1;

va,s = DIVAndrun(mask,(pm,pn),(xi,yi),(x,y),v,(lenx,leny),epsilon2,primal=true)

@test va[2,2] ≈ 1 atol=1e-3


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
