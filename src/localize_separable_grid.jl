# consider to reverse order x,mask,xi
"""
Derive fractional indices on a separable grid.

I = localize_separable_grid(xi,mask,x)

xi and x are a tuples, e.g.
x1,x2 = ndgrid(2 * collect(1:5),collect(1:6))
x = (x1,x2)

Derive fractional indices where xi are the points to localize in the
separable grid `x` (every dimension in independent on other dimension).
The output `I` is an n-by-m array where n number of dimensions and m number of
observations. The correspond element of I is negative if `xi` is outside of
the grid defined by `x`.
"""
function localize_separable_grid(xi,mask,x)

    # n dimension of the problem
    n = length(x)

    # m is the number of arbitrarily distributed observations
    mi = prod(size(xi[1]))

    # sz is the size of the grid
    sz = size(x[1])

    vi = []
    X = []
    I = zeros(n,mi)

    for i=1:n
        push!(X,[x[i][(j-1)*stride(x[i],i) + 1] for j in 1:sz[i]])
        push!(vi,collect(1:sz[i]))
    end

    IJ = ndgrid(vi...)

    for i=1:n
        #@show typeof(X[1])
        #@show (X...)
        #@show size(IJ[i])

        itp = interpolate((X...),IJ[i],Gridded(Linear()))

        # loop over all point
        for j = 1:mi
            I[i,j] = itp[[_[j] for _ in xi]...]
        end
    end

    # handle rounding errors
    # snap to domain bounding box if difference does not exceeds tol
    tol = 50*eps(1.)

    for i=1:n
        # upper bound
        ind = sz[i] .< I[i,:] .<= sz[i] + tol
        I[i,ind] = sz[i]

        # lower bound
        ind = 1 .< I[i,:] .<= 1 + tol
        I[i,ind] = 1
    end

    return I
end


# LocalWords:  indices sz tol

# Copyright (C) 2014, 2016, 2017 Alexander Barth <a.barth@ulg.ac.be>
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
