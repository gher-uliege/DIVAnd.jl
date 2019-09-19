"""
Computes an N-dimensional volume integral

DIVAnd_integral(mask,pmn,fieldin)


# Input:
*  `mask`: mask as usual
*  `pmn` : tuple of metrics as usual
*  `fieldin`: field of the same dimensions as mask and which is integrated over the domain

# Output:
*  `integratedfield`: The integral


"""
function DIVAnd_integral(mask, pmn, fieldin)
    # Dimensions but not checked for consistency
    NDIM = ndims(mask)
    dim = size(mask)
    # Utility array holding volume based on metrics
    volume = zeros(Float64, dim)
    volume[mask] .= 1.0
    for i = 1:NDIM
        volume .= volume ./ pmn[i]
    end

    # Avoiding NaN problems and so on putting field not on the grid to zero
    field = fieldin
    field[.!mask] .= 0.



    #@show volume,sum(volume)
    # Return simple integral estimate
    integratedfield = sum(field .* volume)
    return integratedfield
end

# Copyright (C) 2008-2019 Alexander Barth <barth.alexander@gmail.com>
#                         Jean-Marie Beckers   <JM.Beckers@ulg.ac.be>
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


