# Testing DIVAnd in 2 dimensions with independent verification.

if VERSION >= v"0.7.0-beta.0"
    using Test
else
    using Base.Test
end
using DIVAnd


# grid of background field
mask,(pm,pn),(xi,yi) = DIVAnd_squaredom(2,range(0,stop=1,length=30))

epsilon = 1e-10;

# grid of observations
x,y = ndgrid(range(epsilon,stop = 1-epsilon, length = 20),
             range(epsilon,stop = 1-epsilon, length = 20))
x = x[:]
y = y[:]
v = sin.(x*6) .* cos.(y*6)


lenx = .15;
leny = .15;

epsilon2 = 0.05;

#,err,s
va,s = DIVAndrun(mask,(pm,pn),(xi,yi),(x,y),v,(lenx,leny),epsilon2,primal=true)

#Z = randn(size(s.H,1),size(s.H,1));
Z = eye(size(s.H,1));

ZtHKZ = Z' * (s.H*(s.P * (s.H'* (s.R \ Z))));
WW=s.P * (s.H'* (s.R \ Z)); ZtHKZ2 =  Z'*s.H*WW;

@time ZtHKZ = Z'*s.H*(s.P * (s.H'* (s.R \ Z)));

@time begin
    WW=s.P * (s.H'* (s.R \ Z));
    ZtHKZ2 =  Z'*s.H*WW;
end





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
