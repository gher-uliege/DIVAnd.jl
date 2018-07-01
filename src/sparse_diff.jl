function _sparse_wsum(sz1,m,v1,v2,cyclic = false)

    n1 = prod(sz1)

    sz2 =
        if cyclic
            sz1
        else
            ntuple(i -> (i == m ? sz1[i]-1 : sz1[i] ),length(sz1))
        end

    n2 = prod(sz2)
    n = length(sz1)

    L2 = zeros(Int,n2)
    L2o = zeros(Int,n2)
    L1 = 1:n2
    v = ones(size(L1))

    tmp = zeros(Int,n)

    for i = 1:n2
        lin2sub!(sz2,i,tmp)
        L2[i] = sub2lin(sz1,tmp)

        tmp[m] += 1
        if cyclic
            tmp[m] = mod(tmp[m] - 1,sz1[m]) + 1
        end

        L2o[i] = sub2lin(sz1,tmp)
    end


    S = sparse(
               [L1;     L1;  ],
               [L2;     L2o; ],
               [fill(v1,n2);  fill(v2,n2)], n2 , n1 )

    return S

end


"""
    diffx = sparse_diff(sz1,m,cyclic)

Sparse operator for differentiation along dimension `m` for "collapsed" matrix
of the size `sz1`. `cyclic` is true if domain is cyclic along dimension m.
`false` is the default value
"""
sparse_diff(sz1,m,cyclic = false) = _sparse_wsum(sz1,m,-1.,1.,cyclic)


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
