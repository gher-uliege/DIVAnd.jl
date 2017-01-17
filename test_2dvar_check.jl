# Testing divand in 2 dimensions with independent verification.

using Base.Test
#using divand

# grid of background field
xi,yi = ndgrid(linspace(0,1,10),linspace(0,1,10))

mask = trues(size(xi))
pm = ones(size(xi)) / (xi[2,1]-xi[1,1])
pn = ones(size(xi)) / (yi[1,2]-yi[1,1])

epsilon = 1e-10;

# grid of observations
x,y = ndgrid(linspace(epsilon,1-epsilon,10),linspace(epsilon,1-epsilon,10))
x = x[:]
y = y[:]
v = sin( x*6 ) .* cos( y*6)


lenx = .15;
leny = .15;

lambda = 20;

#,err,s
va,s = divandrun(mask,(pm,pn),(xi,yi),(x,y),v,(lenx,leny),lambda,diagnostics=true,primal=true)

#err = diag(s.P)
iR = inv(full(s.R));
iB = full(s.iB);
H = full(s.H);
sv = s.sv;

iP = iB + H'*iR*H;

P = inv(iP);

xa2 = P* (H'*iR*v[:]);

fi2, = statevector_unpack(sv,xa2);
fi2[~s.mask] = NaN;

@test va ≈ fi2
#@test diag(P) ≈ err[:]


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
