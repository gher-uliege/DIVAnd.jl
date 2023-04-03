# Testing DIVAnd in 2 dimensions with advection.

using Test
using StableRNGs
using Random
using DIVAnd

rng = StableRNG(1234)

# grid of background field
mask, (pm, pn), (xi, yi) = DIVAnd_squaredom(2, range(-1, stop = 1, length = 30))
epsilon2 = 1 / 200
errmean=0

a = 5.0
u = a * yi;
v = -a * xi;


for ND in (10, 1000, 10000)
for len in (0.05, 0.2 ,1.0)
for mymeth in (:auto, :cheap, :precise, :cpme, :scpme, :exact, :aexerr, :diagapp)
for myscale in (true, false)

xdd = rand(rng,ND)
ydd = rand(rng,ND)
fdd = randn(rng,ND)


#@show xdd[1],ND,len,mymeth,myscale,size(mask),size(pm),size(xi),size(yi)


fi11,sl = DIVAndrun(
    mask,
    (pm, pn),
    (xi, yi),
    (xdd, ydd),
    fdd,
    len,
    epsilon2;
    velocity = (u, v),
    alphabc = 0
)



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
	Bscale=myscale,
	rng=StableRNG(123)
	)

	#@show size(errorm),size(fi11)

global errmean=errmean+errorm[10,10]
#@show errmean,errorm[10,10],fi11[10,10],xdd[1],ND,len,mymeth,myscale,methodc

end
end
end
end
@show errmean
@test abs(errmean  - 7.527185124068747) < 0.00000001


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
