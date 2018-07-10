
struct statevector
    mask
    nvar
    numels
    numels_all
    size
    ind
    n
    packed2unpacked
    unpacked2packed
end

function unpack(v,mask)
    tmp = zeros(eltype(v),mask);
    tmp[mask] = v;
    return tmp;
end


"""
Initialize structure for packing and unpacking given their mask.

s = statevector_init((mask1, mask2, ...))

Initialize structure for packing and unpacking
multiple variables given their corresponding land-sea mask.

Input:
  mask1, mask2,...: land-sea mask for variable 1,2,... Sea grid points correspond to one and land grid points to zero.
    Every mask can have a different shape.

Output:
  s: structure to be used with statevector_pack and statevector_unpack.

Note:
see also statevector_pack, statevector_unpack

Author: Alexander Barth, 2009,2017 <a.barth@ulg.ac.be>
License: GPL 2 or later
"""
function statevector_init(masks)

numels = [sum(mask)    for mask in masks]
ind = [0 cumsum(numels)...]

# vector mapping packed indices to unpacked indices
packed2unpacked = [(1:length(mask))[mask] for mask in masks]

# vector mapping unpacked indices packed indices
unpacked2packed = [unpack(1:sum(mask),mask) for mask in masks]

s = statevector(
     [mask for mask in masks],
     length(masks),
     numels,
     [length(mask) for mask in masks],
     [size(mask) for mask in masks],
     ind,
     ind[end],
     packed2unpacked,
     unpacked2packed
)

return s
end

# Copyright (C) 2009,2017 Alexander Barth <a.barth@ulg.ac.be>
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


#  LocalWords:  statevector init GPL
