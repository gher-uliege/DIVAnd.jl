# Testing DIVAnd in 2 dimensions with independent verification.
using DIVAnd

if VERSION >= v"0.7.0-beta.0"
    using Test
else
    using Base.Test
end

# grid of background field
mask,(pm,pn),(xi,yi) = DIVAnd_squaredom(2,range(-1, stop = 1, length = 10))

# make sure that observations are strictly inside the domain
# defined by xi and yi
epsilon = 1e-10;

# grid of observations
x,y = ndgrid(
    range(epsilon,stop = 1-epsilon,length = 10),
    range(epsilon,stop = 1-epsilon,length = 10))

x = x[:]
y = y[:]
v = sin.(6x) .* cos.(6y)

# correlation length
lenx = .15;
leny = .15;

"""naive analysis using full matrices"""
function naive_analysis(s,v)
    iR = inv(Matrix(s.R));
    iB = Matrix(s.iB);
    H = Matrix(s.H);
    iP = iB + H'*iR*H;
    P = inv(iP);
    xa = P* (H'*iR*v[:]);
    return (statevector_unpack(s.sv,xa)[1], statevector_unpack(s.sv,diag(P))[1])
end

# diagonal R with constant diagonal elements

epsilon2 = 0.05;
xa,s = DIVAndrun(mask,(pm,pn),(xi,yi),(x,y),v,(lenx,leny),epsilon2,primal=true)
diagP, = statevector_unpack(s.sv,diag(s.P))
xa_check, diagP_check = naive_analysis(s,v)
@test xa ≈ xa_check
@test diagP ≈ diagP_check

# diagonal R with varying diagonal elements

diagR = sin.(4x) .+ 2
xa,s = DIVAndrun(mask,(pm,pn),(xi,yi),(x,y),v,(lenx,leny),diagR,primal=true)
diagP, = statevector_unpack(s.sv,diag(s.P))
xa_check, diagP_check = naive_analysis(s,v)
@test xa ≈ xa_check
@test diagP ≈ diagP_check


# diagonal R with varying diagonal elements (2)

R = Diagonal(sin.(4x) .+ 2)
xa,s = DIVAndrun(mask,(pm,pn),(xi,yi),(x,y),v,(lenx,leny),R,primal=true)
diagP, = statevector_unpack(s.sv,diag(s.P))
xa_check, diagP_check = naive_analysis(s,v)
@test xa ≈ xa_check
@test diagP ≈ diagP_check

# non-diagonal R
m = length(x)

#R = spdiagm((ones(m-1),4*ones(m),ones(m-1)),(-1,0,1))
R = spdiagm(-1 => ones(m-1), 0 => 4*ones(m), 1 => ones(m-1))

xa,s = DIVAndrun(mask,(pm,pn),(xi,yi),(x,y),v,(lenx,leny),R,primal=true)
diagP, = statevector_unpack(s.sv,diag(s.P))
xa_check, diagP_check = naive_analysis(s,v)
@test xa ≈ xa_check
@test diagP ≈ diagP_check

# Copyright (C) 2014-2017 Alexander Barth <a.barth@ulg.ac.be>
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
