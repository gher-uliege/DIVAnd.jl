# A simple example of DIVAnd in 4 dimensions
# with observations from an analytical function.

using DIVAnd
if VERSION >= v"0.7.0-beta.0"
    using Test
else
    using Base.Test
end
#using PyPlot

# final grid
gridsize = (101,101)

n = length(gridsize)

# observations
xy = ntuple(i -> [0.5],n)
f = [1]


# mask: all points are valid points
# this problem has a simple cartesian metric
# pm is the inverse of the resolution along the 1st dimension,...


mask,pmn,xyi = DIVAnd_rectdom(
    [range(0, stop = 1,length = s) for s in gridsize]...)


sv = statevector((mask,))

# correlation length
#lenxy = ntuple(i -> 1.,n)
lenxy = ntuple(i -> .1,n)
lenxy = (0.1,0.15)

# obs. error variance normalized by the background error variance
epsilon2 = 1.;


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

Ld = [mean(L) for L in lenxy]
nu = ([L.^2 for L in lenxy]...)
D = DIVAnd_laplacian(Val{:sparse},mask,pmn,nu,falses(4))



#1/dx 1/dy ( ∂_x  dy/dx nu  ∂_x  )

zi0 = zeros(gridsize)
zi0[50,50] = 1

ivol = .*(pmn...)

@show sum(zi0 ./ ivol)


zi = zi0[:]



# G(x,x',t) = (4πDt)^(-n/2)  exp( - |x -x'|² / (4Dt))

nmax2 = 10000
α = 0.0001

for niter = 1:nmax2
    zi += α * (D*zi)
end


zi = reshape(zi,gridsize)

@show sum(zi ./ ivol)

xi,yi = xyi



ri = [ sqrt((xi[i,j]-xi[50,50])^2 / Ld[1]^2  + (yi[i,j]-yi[50,50])^2 / Ld[2]^2)  for i=1:gridsize[1], j=1:gridsize[2]]
zit = exp(-ri.^2/(4 * nmax2 * α)) ./ ((4 * π * nmax2 * α)^(n/2) * prod(Ld) * ivol[50,50])

@show maximum(abs(zit-zi))
@test_approx_eq_eps zit zi 1e-4

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
