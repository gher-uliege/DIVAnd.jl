if VERSION >= v"0.7.0-beta.0"
    using Test
else
    using Base.Test
end

mask = rand(10,10) .> .5;
mask_u = rand(9,10) .> .5;
mask_v = rand(10,9) .> .5;

# a couple of point should always be unmasked
mask[3,5] = true
mask_u[4,6] = true
mask_v[5,7] = true


sv = statevector((mask,mask_u,mask_v));
var = rand(10,10);
var[mask.==0] .= 0;

var_u = rand(9,10);
var_u[mask_u.==0] .= 0;

var_v = rand(10,9);
var_v[mask_v.==0] .= 0;

E = pack(sv,(var,var_u,var_v));
Ezeta2,Eu2,Ev2 = unpack(sv,E);

#@show all(Ezeta2 .== var)
#@show maximum(abs(Ezeta2 - var))

@test Ezeta2 ≈ var
@test Eu2 ≈ var_u
@test Ev2 ≈ var_v


data0 = randn(sv.n,5)
Ezeta2,Eu2,Ev2 = unpackens(sv,data0)
data2 = packens(sv,(Ezeta2,Eu2,Ev2))
@test size(Ezeta2,3) == 5
@test data0 ≈ data2


ind = sub2ind(sv,(1,3,5))
@test var[3,5] ≈ E[ind]
@test ind2sub(sv,ind) == (1,3,5)


ind = sub2ind(sv,(2,4,6))
@test var_u[4,6] ≈ E[ind]
@test ind2sub(sv,ind) == (2,4,6)

ind = sub2ind(sv,(3,5,7))
@test var_v[5,7] ≈ E[ind]
@test ind2sub(sv,ind) == (3,5,7)


# mask as Array{Bool,1} instead of BitArray
sv = statevector(([true,false],))
@test sv.n == 1


# Copyright (C) 2009 Alexander Barth <a.barth@ulg.ac.be>
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; If not, see <http://www.gnu.org/licenses/>.


#  LocalWords:  statevector init GPL
