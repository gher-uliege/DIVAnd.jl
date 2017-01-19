using Base.Test

mask = rand(10,10) .> .5;
mask_u = rand(9,10) .> .5;
mask_v = rand(10,9) .> .5;

sv = statevector_init((mask,mask_u,mask_v));
var = rand(10,10);
var[mask.==0] = 0;

var_u = rand(9,10);
var_u[mask_u.==0] = 0;

var_v = rand(10,9);
var_v[mask_v.==0] = 0;

E = statevector_pack(sv,(var,var_u,var_v));
Ezeta2,Eu2,Ev2 = statevector_unpack(sv,E);


#@show all(Ezeta2 .== var)
#@show maximum(abs(Ezeta2 - var))

@test Ezeta2 ≈ var
@test Eu2 ≈ var_u
@test Ev2 ≈ var_v



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
