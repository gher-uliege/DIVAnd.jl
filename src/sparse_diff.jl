"""
Sparse operator for differentiation.

diffx = sparse_diff(sz1,m,cyclic)

Sparse operator for differentiation along dimension m for "collapsed" matrix
of the size sz1.

Input:
  sz1: size of rhs
  m: dimension to differentiate
  cyclic: true if domain is cyclic along dimension m. False is the
  default value
"""

function sparse_diff(sz1,m,cyclic = false)

n1 = prod(sz1)

# make a copy and convert tuple to array
sz2 = collect(sz1)

if !cyclic
    sz2[m] = sz2[m]-1
end

n2 = prod(sz2)
n2i = Int(n2)
n = length(sz1)


strides1 = [1,cumprod(collect(sz1))[1:end-1]...]

L2 = Array{Int64,1}(n2)


# L2 = [i for i = 1:n1 if ind2sub(sz1,i)[m] != sz1[m] ]
# L2o = L2 + strides1[m]


vi = [collect(1:sz2[i]) for i = 1:n]
#IJ = [_[:] for _ in ndgrid(vi...)]

#IJ = [Array{Int64,1}(n2) for i = 1:n]
IJ = [zeros(Int64,Int(n2)) for i = 1:n]

    s = 1
    for i=1:n
        v = vi[i]
        snext = s*sz2[i]
        ndgrid_fill(IJ[i], v, s, snext)
        s = snext
    end



#L1 = [1:Int(n2);]
L1 = [1:n2i;]
L2 = sub2ind(sz1,IJ...)::Array{Int64,1}
one = ones(size(L1))

IJ[m] = IJ[m] + 1

if cyclic
     IJ[m] = mod(IJ[m]-1,sz1[m])+1
end

L2o = sub2ind(sz1,IJ...)::Array{Int64,1}

S = sparse(
    [L1;     L1;  ],
    [L2;     L2o; ],
    [-one;   one  ], n2 , n1 )

return S

end

# Copyright (C) 2009,2012,2017 Alexander Barth <a.barth@ulg.ac.be>
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


# LocalWords:  diffx sz
