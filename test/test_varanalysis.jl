# A simple example of divand in 4 dimensions
# with observations from an analytical function.

using divand
using Base.Test
#using PyPlot

# final grid
gridsize = (101,101)
#gridsize = (101,101,101)

n = length(gridsize)

# observations
xy = ntuple(i -> [0.5],n)
f = [1.]


# mask: all points are valid points
# this problem has a simple cartesian metric
# pm is the inverse of the resolution along the 1st dimension,...


mask,pmn,xyi = divand_rectdom([linspace(0,1,s) for s in gridsize]...)


sv = statevector((mask,))

# correlation length
#lenxy = ntuple(i -> 1.,n)
lenxy = ntuple(i -> .1,n)

# obs. error variance normalized by the background error variance
epsilon2 = 1.;


# tolerance on the gradient A x - b
tol = 1e-1
tol = 1e-5



@time xa,s = varanalysis(mask,pmn,xyi,xy,f,lenxy,epsilon2; tol = tol)
@show maximum(xa) 

@test_approx_eq_eps maximum(xa) 0.5 1e-3

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
