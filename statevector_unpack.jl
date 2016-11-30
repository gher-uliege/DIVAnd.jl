# Unpack a vector into different variables under the control of a mask.
#
# [var1, var2, ...] = statevector_unpack(s,x)
# [var1, var2, ...] = statevector_unpack(s,x,fillvalue)
#
# Unpack the vector x into the different variables var1, var2, ...
#
# Input:
#   s: structure generated by statevector_init.
#   x: vector of the packed elements. The size of this vector is the number of elements equal to 1
#     in all masks.
#
# Optional input parameter:
#   fillvalue: The value to fill in var1, var2,... where the masks correspond to a land grid point. The default is zero.
#
# Output:
#   var1, var2,...: unpacked variables.
#
# Notes:
# If x is a matrix, then the second dimension is assumed 
# to represent the different ensemble members. In this case,
# var1, var2, ... have also an additional trailing dimension.

# Author: Alexander Barth, 2009 <a.barth@ulg.ac.be>
# License: GPL 2 or later

function statevector_unpack(s,x,fillvalue = 0)


k = size(x,2)

out = []

for i=1:s.nvar
  v = zeros(s.numels_all[i],k)
  v[:] = fillvalue
  
  ind = find(s.mask[i])

  v[ind,:] = x[s.ind[i]+1:s.ind[i+1],:]
  
  push!(out,reshape(v,([s.size[i]... k]...)))
end


out
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


# LocalWords:  statevector fillvalue init GPL varargout
