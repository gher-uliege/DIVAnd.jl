# A simple example of DIVAnd in 4 dimensions
# with observations from an analytical function.

using DIVAnd
using Test
using LinearAlgebra

# final grid
gridsize = (101, 101)

n = length(gridsize)

# observations
xy = ntuple(i -> [0.], n)
f = [2.]


# mask: all points are valid points
# this problem has a simple cartesian metric
# pm is the inverse of the resolution along the 1st dimension,...

mask, pmn, xyi = DIVAnd_rectdom([range(-1, stop = 1, length = s) for s in gridsize]...)


sv = statevector((mask,))

# correlation length
len = ntuple(i -> 0.2, n)

# obs. error variance normalized by the background error variance
epsilon2 = 1.;


alpha = [1, 2, 1]
fi, s = DIVAndrun(mask, pmn, xyi, xy, f, len, epsilon2, alpha = alpha)
mu, K, len_scale = DIVAnd_kernel(n, alpha);
# xy is a tuple with the coordinates in every dimensions
fit = [K(len_scale * norm([xy...] ./ [len...])) for xy in zip(xyi...)]
@test fi ≈ fit rtol = 1e-2


alpha = [1, 3, 3, 1]
fi, s = DIVAndrun(mask, pmn, xyi, xy, f, len, epsilon2, alpha = alpha)
mu, K, len_scale = DIVAnd_kernel(n, alpha);
# xy is a tuple with the coordinates in every dimensions
fit = [K(len_scale * norm([xy...] ./ [len...])) for xy in zip(xyi...)]
@test fi ≈ fit rtol = 1e-2


alpha = [1, 0, 1]
fi, s = DIVAndrun(mask, pmn, xyi, xy, f, len, epsilon2, alpha = alpha, scale_len = false)
mu, K, len_scale = DIVAnd_kernel(n, alpha);
# xy is a tuple with the coordinates in every dimensions
fit = [K(len_scale * norm([xy...] ./ [len...])) for xy in zip(xyi...)]
@test fi ≈ fit rtol = 0.5




# Copyright (C) 2014, 2017 Alexander Barth <a.barth@ulg.ac.be>
#                          Jean-Marie Beckers <JM.Beckers@ulg.ac.be>
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
