# Initialize structure for packing and unpacking given their mask.
#
# s = statevector_init(mask1, mask2, ...)
#
# Initialize structure for packing and unpacking
# multiple variables given their corresponding land-sea mask.
#
# Input:
#   mask1, mask2,...: land-sea mask for variable 1,2,... Sea grid points correspond to one and land grid points to zero.
#     Every mask can have a different shape.
#
# Output:
#   s: structure to be used with statevector_pack and statevector_unpack.
#
# Note:
# see also statevector_pack, statevector_unpack

# Author: Alexander Barth, 2009 <a.barth@ulg.ac.be>
# License: GPL 2 or later

type statevector
    mask
    nvar
    numels
    numels_all
    size
    ind
    n
end

function statevector_init(masks)

numels = [sum(_)    for _ in masks]
ind = [0 cumsum(numels)...]

s = statevector(
     [_ for _ in masks],
     length(masks),
     numels,
     [length(_) for _ in masks],
     [size(_) for _ in masks],
     ind,
     ind[end]
     )



# s.nvar = nargin;


# for i=1:nargin
#   mask = varargin{i};
#   s.mask{i} = mask;
#   s.numels(i) = sum(mask(:) == 1);
#   s.numels_all(i) = numel(mask);
#   s.size{i} = size(mask);
# end

# s.ind = [0 cumsum(s.numels)];

# s.n = s.ind(end);

s
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


#  LocalWords:  statevector init GPL
