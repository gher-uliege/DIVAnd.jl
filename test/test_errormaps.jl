# Testing DIVAnd in 2 dimensions with advection.
using DIVAnd
using Test
using StableRNGs
using Random

# grid of background field
mask, (pm, pn), (xi, yi) = DIVAnd_squaredom(2, range(-1, stop = 1, length = 30))
epsilon2 = 1 / 200
errmean=0

a = 5;
u = a * yi;
v = -a * xi;


for ND in (10, 1000, 10000)
for len in (0.05, 0.2 ,1.0)
for mymeth in (:auto, :cheap, :precise, :cpme, :scpme, :exact, :aexerr, :diagapp)
for myscale in (true, false)

xdd = rand(ND)
ydd = rand(ND)
fdd = randn(ND)





fil, sl = DIVAndrun(
    mask,
    (pm, pn),
    (xi, yi),
    (xdd, ydd),
    fdd,
    len,
    epsilon2;
    velocity = (u, v),
    alphabc = 0,
);

errorm,methodc=DIVAnd_errormap(
    mask,
    (pm, pn),
    (xi, yi),
    (xdd, ydd),
    fdd,
    len,
    epsilon2,sl;
    velocity = (u, v),
    alphabc = 0,
	method=mymeth,
	Bscale=myscale
	)
	
global errmean=errmean+errorm[10,10]

end
end
end
end

@test errmean  ≈ 6.951700275663623


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
