# Sparse operator for trimming.
#
# T = sparse_trim(sz1,m)
#
# Create a sparse operator which trim first and last row (or column) in 
# The field is a "collapsed" matrix of the size sz1.
#
# Input:
#   sz1: size of rhs
#   m: dimension to trim

function sparse_trim(sz1,m)

n1 = prod(sz1)
sz2 = collect(sz1)
sz2[m] = sz2[m]-2
n2 = prod(sz2)
n = length(sz1)

# L1

vi = [collect(1:sz2[i]) for i = 1:n]
IJ = [_[:] for _ in ndgrid(vi...)]

L1 = sub2ind((sz2...),IJ...)

IJ[m]=IJ[m]+1

L2 = sub2ind(sz1,IJ...)

one = ones(size(L1))

sparse(L1, L2, one, n2, n1)

end

# Copyright (C) 2009,2016 Alexander Barth <a.barth@ulg.ac.be>
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

