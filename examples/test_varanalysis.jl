# A simple example of divand in 4 dimensions
# with observations from an analytical function.

using divand
#using PyPlot

# final grid
gridsize = (101,101)

ndims = length(gridsize)

# observations
xy = ntuple(i -> [0.5],ndims)
f = [1]


# mask: all points are valid points
# this problem has a simple cartesian metric
# pm is the inverse of the resolution along the 1st dimension,...


mask,pmn,xyi = divand_rectdom([linspace(0,1,s) for s in gridsize]...)


sv = statevector((mask,))

# correlation length
lenxy = ntuple(i -> 0.2,ndims)

# obs. error variance normalized by the background error variance
epsilon2 = 1;


# tolerance on the gradient A x - b
tol = 1e-1
tol = 1e-5

# http://www.rpgroup.caltech.edu/~natsirt/aph162/diffusion_old.pdf
# equation 60
# n-dimensional Green’s function
# ∂c  
# --  =  ∇ ⋅ (D ∇ c)
# ∂t 

# G(x,x',t) = (4πDt)^(-n/2)  exp( - |x -x'|² / (4Dt))


@time xa = varanalysis(mask,pmn,xyi,xy,f,lenxy,epsilon2; tol = tol)


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
