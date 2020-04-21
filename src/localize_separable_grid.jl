
# For floating point coordinates, the tolerance is 50 ϵ,
# for integer and other types it is 0.
tolerance(::Type{T}) where T <: AbstractFloat = 50*eps(T)
tolerance(T) = 0

function _localize_separable_grid_cyclic!(
    xi::NTuple{n,AbstractArray},
    x::NTuple{n,AbstractArray{T,n}},
    iscyclic,
    I,
) where {n} where {T}

    # m is the number of arbitrarily distributed observations
    mi = prod(size(xi[1]))

    # sz is the size of the grid
    sz = size(x[1])

    # extent array
    vi = ntuple(i -> 1:(sz[i]+iscyclic[i]), Val(n))
    IJ = ndgrid(vi...)::NTuple{n,Array{Int,n}}

    X = ntuple(i -> Vector{T}(undef, sz[i] + iscyclic[i]), Val(n))
    for i = 1:n
        for j = 1:sz[i]
            X[i][j] = x[i][(j-1)*stride(x[i], i)+1]
        end

        if iscyclic[i]
            X[i][end] = 2 * X[i][end-1] - X[i][end-2]
        end
    end

    scheme = ntuple(i -> (iscyclic[i] ? Periodic() : Throw()), Val(n))

    for i = 1:n
        itp = extrapolate(interpolate(X, IJ[i], Gridded(Linear())), scheme)

        # loop over all point
        for j = 1:mi
            try
                xind = NTuple{n,T}(getindex.(xi, j))
                I[i, j] = itp(xind...)
            catch
                # extrapolation
                I[i, j] = -1
            end
        end
    end

    return I
end

# consider to reverse order x,mask,xi
"""
Derive fractional indices on a separable grid.

    I = localize_separable_grid(xi,mask,x)

xi is a tuple of vectors and x and tuple of n-dimensional arrays, e.g.

x1,x2 = ndgrid(2 * collect(1:5),collect(1:6))
x = (x1,x2)

Derive fractional indices where xi are the points (typical discrete observations) to localize in the
separable grid `x` (every dimension in independent on other dimension).
The output `I` is an n-by-m array where n number of dimensions and m number of
observations. The corresponding element of I is negative if `xi` is outside of
the grid defined by `x`.
"""
function localize_separable_grid(
    xi::NTuple{n,AbstractArray},
    mask::AbstractArray{Bool,n},
    x::NTuple{n,AbstractArray{T,n}},
    iscyclic = falses(n),
) where {n} where {T}

    # m is the number of arbitrarily distributed observations
    mi = prod(size(xi[1]))

    # sz is the size of the grid
    sz = size(x[1])

    I = zeros(n, mi)
    X = ntuple(i -> Vector{T}(undef, sz[i]), Val(n))
    vi = ntuple(i -> 1:sz[i], Val(n))

    for i = 1:n
        for j = 1:sz[i]
            X[i][j] = x[i][(j-1)*stride(x[i], i)+1]
        end
    end

    IJ = ndgrid(vi...)::NTuple{n,Array{Int,n}}

    if any(iscyclic)
        _localize_separable_grid_cyclic!(xi, x, iscyclic, I)
    else
        for i = 1:n
            # https://github.com/JuliaMath/Interpolations.jl/issues/237
            itp = extrapolate(interpolate(X, IJ[i], Gridded(Linear())), -1.0)

            # loop over all point
            for j = 1:mi
                xind = NTuple{n,T}(getindex.(xi, j))
                I[i, j] = itp(xind...)
            end
        end
    end

    indices = I
    #JLD2.@save "/tmp/loc.jld2" X indices IJ xi x
    #@show extrema(I),size(I)
    #@show I[:,1]

    # handle rounding errors
    # snap to domain bounding box if difference does not exceeds tol
    tol = tolerance(T)

    for i = 1:n
        # upper bound
        ind = sz[i] .< I[i, :] .<= sz[i] + tol
        I[i, ind] .= sz[i]

        # lower bound
        ind = 1 .< I[i, :] .<= 1 + tol
        I[i, ind] .= 1
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
