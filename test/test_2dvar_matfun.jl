# Testing DIVAnd in 2 dimensions

if VERSION >= v"0.7.0-beta.0"
    using Test
else
    using Base.Test
end

# number of dimensions
n = 2

# function to interpolate
fun(xy...) = .*([cos.(6*x) for x in xy]...)

# grid of background field
mask,pmn,xyi = DIVAnd_squaredom(2,range(0, stop = 1, length = 20))

epsilon = 1e-10;

# grid of observations
xy = ndgrid([range(epsilon,stop = 1-epsilon,length = 10) for i = 1:n]...)
v = fun([x[:] for x in xy]...);

# correlation length
len = .15;

epsilon2 = 1.

# tolerance on the gradient A x - b
tol = 1e-4
# tolerance on the result x
tolres = 1e-3

kwargs = [(:tol, tol),(:maxit,10000)]
# type of operators Val{:sparse} or Val{:MatFun}
operatortype = Val{:sparse}
#operatortype = Val{:MatFun}

# iterative (without preconditioner)
vas,s_np = DIVAndrun(mask,pmn,xyi,xy,v,len,epsilon2;
                     kwargs..., inversion=:pcg,operatortype=Val{:sparse})

#@show typeof(s_np.iB)

vamf,s_np = DIVAndrun(mask,pmn,xyi,xy,v,len,epsilon2;
                      kwargs..., inversion=:pcg,operatortype=Val{:MatFun})

#@show typeof(s_np.iB)

#@show s_np.niter

@test vas ≈ vamf atol=2e-4


# Copyright (C) 2017 Alexander Barth <a.barth@ulg.ac.be>
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
