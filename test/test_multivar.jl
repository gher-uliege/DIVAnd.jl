# A simple example of DIVAnd in 2 dimensions
# with observations from an analytical function.

using Test
using StableRNGs
using Random
using DIVAnd

# true error variance of observation
epsilon2_true = 1.0

rng = StableRNG(1234)
#rng = MersenneTwister(1234)

# observations
nobs = 99
x = rand(rng,nobs);
y = rand(rng,nobs);
v = rand(rng,[1,2],nobs);

# true length-scale
len_true = 0.5
f = (1.5 .-v ).*sin.(π * x / len_true) .* cos.(π * y / len_true);

# add noise
f = f + sqrt(epsilon2_true) * randn(rng,nobs);

# final grid
mask, (pm, pn,pv), (xi, yi,vi) =
    DIVAnd_rectdom(range(0, stop = 1, length = 14), range(0, stop = 1, length = 13),1:2)

# correlation length (first guess)
len = 0.1;

# obs. error variance normalized by the background error variance (first guess)
epsilon2 = 2.0

# loop over all methods

fi,s=DIVAnd_multivarEOF(mask,(pm,pn,pv),(xi,yi,vi),(x,y,v),f,len,epsilon2;velocity=(ones(size(mask)),ones(size(mask)),zeros(size(mask))))


@test fi[10,10,2] ≈  0.15748780548566835

fi,s=DIVAnd_multivarEOF(mask,(pm,pn,pv),(xi,yi,vi),(x,y,v),f,(len,len,0.),epsilon2;velocity=(ones(size(mask)),ones(size(mask)),zeros(size(mask))),eof=[-1.,1.])

@test fi[10,10,2] ≈ 0.1177009988865455

fi,s=DIVAnd_multivarJAC(mask,(pm,pn,pv),(xi,yi,vi),(x,y,v),f,len,epsilon2;velocity=(ones(size(mask)),ones(size(mask)),zeros(size(mask))))


@test fi[10,10,2] ≈ 0.051194689308468086

fi,s=DIVAnd_multivarJAC(mask,(pm,pn,pv),(xi,yi,vi),(x,y,v),f,(len,len,0.),epsilon2;velocity=(ones(size(mask)),ones(size(mask)),zeros(size(mask))))

@test fi[10,10,2] ≈ 0.051194689308468086


# Copyright (C) 2014, 2017, 2018
#   Alexander Barth <a.barth@ulg.ac.be>
#   Jean-Marie Beckers <JM.Beckers@ulg.ac.be>
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
