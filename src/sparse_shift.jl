# Sparse operator shifting a field in a given dimension.
#
# function S = sparse_shift(sz1,m,cyclic)
#
# Sparse operator shifting a field in the dimension m. The field is a
# "collapsed" matrix of the size sz1.
#
# Input:
#   sz1: size of rhs
#   m: dimension to shift
#   cyclic: true if domain is cyclic along dimension m. False is the
#     default value

function sparse_shift(sz1,m,cyclic = false)

    n1 = prod(sz1)
    sz2 = collect(sz1)

    if !cyclic
        sz2[m] = sz2[m]-1
    end

    n2 = prod(sz2)
    n = length(sz1)

    vi = [collect(1:sz2[i]) for i = 1:n]
    IJ = [vii[:] for vii in ndgrid(vi...)]

    L1 = 1:n2
    one = ones(size(L1))

    IJ[m] = IJ[m] + 1

    if cyclic
        IJ[m] = mod.(IJ[m]-1,sz1[m])+1
    end

    L2 = sub2ind(sz1,IJ...)
    S = sparse(L1,L2,one,n2,n1)

    return S

end

# Copyright (C) 2012-2017 Alexander Barth <a.barth@ulg.ac.be>
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; If not, see <http://www.gnu.org/licenses/>.
