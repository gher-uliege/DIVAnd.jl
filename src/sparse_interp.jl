"""
    H,out = sparse_interp(mask,I)

Create interpolation matrix from `mask` and fractional indexes `I`.

Input:
  mask: 0 invalid and 1 valid points (n-dimensional array)
  I: fractional indexes (2-dim array n by mi, where mi is the number of points to interpolate)
Ouput:
  H: sparse matrix with interpolation coefficients
  out: true if value outside of grid
  outbbox: 1 if outise bouding box
  onland: 1 if point touches land (where mask == 0)
"""
function sparse_interp(mask,I,iscyclic = falses(size(I,1)))
    if (ndims(mask) != size(I,1)) || (ndims(mask) != length(iscyclic))
        error("sparse_interp: inconsistent arguments")
    end

    sz = size(mask)

    # n dimension of the problem
    n = size(I,1)
    m = prod(sz)

    # mi is the number of arbitrarly distributed observations
    mi = size(I,2)

    # handle cyclic dimensions
    @inbounds for j = 1:mi
        for i = 1:n
            if iscyclic[i]
                # bring I(i,:) inside the interval [1 sz(i)+1[
                # since the i-th dimension is cyclic

                I[i,j] = mod(I[i,j]-1,sz[i])+1
            end
        end
    end

    # integer index
    ind = floor.(Int,I)
    out = falses(mi)

    @inbounds for j = 1:mi
        for i = 1:n
            if !iscyclic[i]
                # make a range check only for non-cyclic dimension
                out[j] = out[j] || !(1 <= I[i,j] <= sz[i])

                # handle border cases
                if I[i,j] == sz[i]
                    ind[i,j] = sz[i]-1
                end
            end
        end
    end

    outbbox = copy(out)

    # number of obs. inside
    mip = mi - sum(out)

    si = zeros(Int,2^n, mip)
    sj = ones(Int,2^n, mip)
    ss = ones(2^n, mip)
    k = 0

    @inbounds for k2 = 1:mi
        # consider now only points inside
        if !out[k2]
            k=k+1

            # loop over all corner of hypercube
            for i=1:2^n
                si[i,k] = k2

                # stride for array
                scale = 1

                # loop over all dimensions
                for j=1:n
                    bit = (i >> (j-1)) & 1

                    # the j-index of the i-th corner has the index ip
                    # this index ip is zero-based
                    ip = ind[j,k2] + bit - 1

                    # ip must be [0 and sz[j]-1] (zero-based)
                    # we know aleady that the point is inside the domain
                    # so, if it is outside this range then it is because of periodicity
                    ip = mod(ip,sz[j])

                    sj[i,k] = sj[i,k] + scale * ip

                    # interpolation coefficient
                    α = I[j,k2] - ind[j,k2]

                    if bit == 0
                        ss[i,k] = ss[i,k] * (1-α)
                    else
                        ss[i,k] = ss[i,k] * α
                    end

                    scale = scale * size(mask,j)
                end
            end

            inside = false
            sum_ss = zero(eltype(ss))

            # point is inside if there is any sea surouding grid point (with a
            # weight different from zero)
            @inbounds for i=1:2^n
                if mask[sj[i,k]] && (ss[i,k] != 0)
                    sum_ss += ss[i,k]
                    inside = true
                end
            end

            # inside is only true if sum_ss is different from zero
            if inside
                @inbounds for i=1:2^n
                    # interpolation weight are scaled such that their sum is 1
                    # but taking only sea points into account
                    if mask[sj[i,k]]
                        ss[i,k] = ss[i,k]/sum_ss
                    else
                        ss[i,k] = 0
                    end
                end
            end

            out[k2] = !inside

        end
    end

    @debug "number of point outside $(sum(out)) or  $(100 * sum(out) / length(out)) % "
    H = sparse(si[:],sj[:],ss[:],mi,m)
    return H,out,outbbox

end

"""
sparse_interp(x,mask,xi)
Interpolate from x onto xi
"""
function sparse_interp_g(x,mask,xi)
    I = localize_separable_grid(xi,mask,x)
    return sparse_interp(mask,I)
end


# LocalWords:  interp

# Copyright (C) 2004, 2016 Alexander Barth <a.barth@ulg.ac.be>
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
