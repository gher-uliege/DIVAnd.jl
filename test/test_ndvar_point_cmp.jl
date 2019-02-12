# A simple example of DIVAnd in 4 dimensions
# with observations from an analytical function.

using DIVAnd
if VERSION >= v"0.7.0-beta.0"
    using Test
else
    using Base.Test
end

# final grid
gridsize = (101,101)

n = length(gridsize)

# observations
xy = ntuple(i -> [0.],n)
f = [2.]


# mask: all points are valid points
# this problem has a simple cartesian metric
# pm is the inverse of the resolution along the 1st dimension,...

mask,pmn,xyi = DIVAnd_rectdom([range(-1,stop=1,length=s) for s in gridsize]...)


sv = statevector((mask,))

# correlation length
len = ntuple(i -> 0.2,n)

# obs. error variance normalized by the background error variance
epsilon2 = 1.;


alpha = [1,2,1]
fi,s = DIVAndrun(mask,pmn,xyi,xy,f,len,epsilon2,alpha = alpha)


coeff_laplacian = zeros(ndims(mask))
coeff_derivative2 = ones(ndims(mask))

fi2,s2 = DIVAndrun(mask,pmn,xyi,xy,f,len,epsilon2,alpha = alpha,
                   coeff_laplacian = coeff_laplacian,
                   coeff_derivative2 = coeff_derivative2
                   )


@test fi ≈ fi2 rtol=0.3



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
