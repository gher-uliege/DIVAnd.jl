
oper_diff(::Type{Val{:sparse}},sz1,m,cyclic = false) = sparse_diff(sz1,m,cyclic)
oper_diff(::Type{Val{:MatFun}},sz1,m,cyclic = false) = matfun_diff(sz1,m,cyclic)

"""
Sparse operator for differentiation.

diffx = matfun_diff(sz1,m,cyclic)

Sparse operator for differentiation along dimension m for "collapsed" matrix
of the size sz1.

Input:
  sz1: size of rhs
  m: dimension to differentiate
  cyclic: true if domain is cyclic along dimension m. False is the
  default value
"""

function matfun_diff(sz1,m,cyclic = false)

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

return MatFun(S)

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
# Sparse operator shifting a field in a given dimension.
#
# function S = matfun_shift(sz1,m,cyclic)
#
# Sparse operator shifting a field in the dimension m. The field is a
# "collapsed" matrix of the size sz1.
#
# Input:
#   sz1: size of rhs
#   m: dimension to shift
#   cyclic: true if domain is cyclic along dimension m. False is the
#     default value

function matfun_shift(sz1,m,cyclic = false)

n1 = prod(sz1)
sz2 = collect(sz1)

if !cyclic
    sz2[m] = sz2[m]-1
end

n2 = prod(sz2)
n = length(sz1)

vi = [collect(1:sz2[i]) for i = 1:n]
IJ = [_[:] for _ in ndgrid(vi...)]

L1 = 1:n2
one = ones(size(L1))

IJ[m] = IJ[m] + 1

if cyclic
    IJ[m] = mod(IJ[m]-1,sz1[m])+1
end

L2 = sub2ind(sz1,IJ...)
S = sparse(L1,L2,one,n2,n1)

return MatFun(S)

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
# Sparse operator for staggering.
#
# S = matfun_stagger(sz1,m,cyclic)
#
# Create a sparse operator for staggering a field in dimension m.
# The field is a "collapsed" matrix of the size sz1.
#
# Input:
#   sz1: size of rhs
#   m: dimension to stagger
#   cyclic: true if domain is cyclic along dimension m. False is the
#   default value

#function matfun_stagger(sz1,m,cyclic = false)::SparseMatrixCSC{Float64,Int64}
function matfun_stagger(sz1,m,cyclic = false)

n1 = prod(sz1)

sz2 = collect(sz1)

if !cyclic
    sz2[m] = sz2[m]-1
end
n2 = prod(sz2)

n = length(sz1)

vi = []

for i=1:n
    push!(vi,collect(1:sz2[i]))
end


tmp = ndgrid(vi...)
IJ = []

for i=1:n
    push!(IJ,tmp[i][:])
end

L1 = 1:n2
L2 = sub2ind(sz1,IJ...)
v = ones(size(L1))/2

IJ[m] = IJ[m] + 1

if cyclic
    IJ[m] = mod(IJ[m]-1,sz1[m])+1
end

L2o = sub2ind(sz1,IJ...)

S = sparse(
       [L1;  L1;  ],
       [L2;  L2o; ],
       [v;   v    ], n2 , n1 );


return MatFun(S)

end

# Copyright (C) 2009 Alexander Barth <a.barth@ulg.ac.be>
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
# Sparse operator for trimming.
#
# T = matfun_trim(sz1,m)
#
# Create a sparse operator which trim first and last row (or column) in
# The field is a "collapsed" matrix of the size sz1.
#
# Input:
#   sz1: size of rhs
#   m: dimension to trim

#function matfun_trim(sz1,m)::SparseMatrixCSC{Float64,Int64}
function matfun_trim(sz1,m)

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

S = sparse(L1, L2, one, n2, n1)

@show typeof(S)
return MatFun(S)

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
