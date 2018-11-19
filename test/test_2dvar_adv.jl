# Testing DIVAnd in 2 dimensions with advection.
using DIVAnd
if VERSION >= v"0.7.0-beta.0"
    using Test
else
    using Base.Test
end

# grid of background field
mask,(pm,pn),(xi,yi) = DIVAnd_squaredom(2,range(-1, stop = 1, length = 30))

x = [.4]
y = [.4]
f = [1.]

a = 5;
u = a*yi;
v = -a*xi;
epsilon2 = 1/200
len = 0.2

fi,s = DIVAndrun(mask,(pm,pn),(xi,yi),(x,y),f,len,epsilon2;velocity = (u,v),alphabc=0);

@test abs(fi[18,24] - 0.8993529043140029) < 1e-2



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
