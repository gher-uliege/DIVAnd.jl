

oper_pack(::Type{Val{:sparse}},mask) = sparse_pack(mask)
oper_pack(::Type{Val{:MatFun}},mask) = matfun_pack(mask)

oper_diag(::Type{Val{:sparse}},d) = sparse_diag(d)
oper_diag(::Type{Val{:MatFun}},d) = matfun_diag(d)

oper_trim(::Type{Val{:sparse}},sz1,m) = sparse_trim(sz1,m)
oper_trim(::Type{Val{:MatFun}},sz1,m) = matfun_trim(sz1,m)

for fun in [:diff,:shift,:stagger]
    @eval begin
        $(Symbol("oper_" * string(fun)))(::Type{Val{:sparse}},sz1,m,cyclic = false) = $(Symbol("sparse_" * string(fun)))(sz1,m,cyclic)
        $(Symbol("oper_" * string(fun)))(::Type{Val{:MatFun}},sz1,m,cyclic = false) = $(Symbol("matfun_" * string(fun)))(sz1,m,cyclic)
    end
end

matfun_diag(d) = MatFun((size(d,1),size(d,1)), x -> d.*x, x -> d.*x)

function matfun_pack(mask)
    MatFun(sparse_pack(mask))
end

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

# sz2 size of the resulting array
sz2 = ntuple(i -> (i == m && !cyclic ? sz1[i]-1 : sz1[i]), length(sz1))

function fun(x)
    x = reshape(x,sz1)

    if !cyclic
        ind1 = [ (i == m ? (2:sz1[i]) : (1:sz1[i])) for i = 1:length(sz1)]
        ind2 = [ (i == m ? (1:sz1[i]-1) : (1:sz1[i])) for i = 1:length(sz1)]
        return (x[ind1...] - x[ind2...])[:]
    else
        ind = [ (i == m ? [2:sz1[i]; 1] : (1:sz1[i])) for i = 1:length(sz1)]
        return (x[ind...] - x)[:]
    end
end

# adjoint
function funt(x)
    x = reshape(x,sz2)
    @show "here"

    if !cyclic
        ind0 = [ (i == m ? 1            : (1:sz2[i])) for i = 1:length(sz2)]
        ind1 = [ (i == m ? (2:sz2[i])   : (1:sz2[i])) for i = 1:length(sz2)]
        ind2 = [ (i == m ? (1:sz2[i]-1) : (1:sz2[i])) for i = 1:length(sz2)]
        ind3 = [ (i == m ? sz2[i]       : (1:sz2[i])) for i = 1:length(sz2)]
        #return (x[ind1...] - x[ind2...])[:]
        @show "here 1D"
        return cat(m,-x[ind0...],x[ind2...]-x[ind1...],x[ind3...])[:]
    else
        ind = [ (i == m ? [sz1[i]; 1:sz1[i]-1] : (1:sz1[i])) for i = 1:length(sz1)]
        return (x[ind...] - x)[:]
    end
end
return MatFun((prod(sz2),prod(sz1)),fun,funt)

end


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

# sz2 size of the resulting array
sz2 = ntuple(i -> (i == m && !cyclic ? sz1[i]-1 : sz1[i]), length(sz1))


function fun(x)
    x = reshape(x,sz1)

    if !cyclic
        ind = [ (i == m ? (2:sz1[i]) : (1:sz1[i])) for i = 1:length(sz1)]
        return x[ind...][:]
    else
        ind = [ (i == m ? [2:sz1[i]; 1] : (1:sz1[i])) for i = 1:length(sz1)]
        return x[ind...][:]
    end
end

# adjoint
function funt(x)
    x = reshape(x,sz2)

    if !cyclic
        @show "here shift nc"
        sz0 = ([ (i == m ? 1 : sz2[i]) for i = 1:length(sz2)]...)
        return cat(m,zeros(eltype(x),sz0),x)[:]
    else
        @show "here shift c"
        ind = [ (i == m ? [sz1[i]; 1:sz1[i]-1] : (1:sz1[i])) for i = 1:length(sz1)]
        return x[ind...][:]
    end
end

return MatFun((prod(sz2),prod(sz1)),fun,funt)



#return MatFun(S)

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
