# Testing DIVAnd in 2 dimensions with data outside of the domain

#using DIVAnd

# grid of background field
mask, (pm, pn), (xi, yi) = DIVAnd_squaredom(2, range(0, stop = 1, length = 100))

epsilon = 1e-10;

# grid of observations
x = [0.5, 2.2]
y = [0.5, .5]
v = [1., 1.]


lenx = .15;
leny = .15;

epsilon2 = 0.05;

# analysis with all values
va, s = DIVAndrun(mask, (pm, pn), (xi, yi), (x, y), v, (lenx, leny), epsilon2)

# analysis with only values inside the domain
inside = 1:1
va2, s = DIVAndrun(
    mask,
    (pm, pn),
    (xi, yi),
    (x[inside], y[inside]),
    v[inside],
    (lenx, leny),
    epsilon2,
)

@test va == va2


x = [0.5, .9]
y = [0.5, .5]
v = [1., 1.]
mask[end-40:end, :] .= false
va, s = DIVAndrun(mask, (pm, pn), (xi, yi), (x, y), v, (lenx, leny), epsilon2)

inside = 1:1
va2, s = DIVAndrun(
    mask,
    (pm, pn),
    (xi, yi),
    (x[inside], y[inside]),
    v[inside],
    (lenx, leny),
    epsilon2,
)

@test va[mask] == va2[mask]


# Copyright (C) 2017 Alexander Barth <a.barth@ulg.ac.be>
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
