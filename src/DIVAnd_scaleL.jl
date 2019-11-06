"""
Computes a relative length based on the mask, metrics and a density field, typically measuring the observation density calculated with DIVAnd_heatmap

lambda = DIVAnd_scaleL(mask,pmn,dens)


# Input:
*  `mask`: mask as usual
*  `pmn` : tuple of metrics as usual
*  `dens`: field of the same dimensions as mask. Higher values of dens will result in lower values of lambda.

# Output:
*  `lambda`: field to be applied to a reference length field. Values are around 1 so some regions will have smaller L and some higher L


"""
function DIVAnd_scaleL(mask, pmn, dens)
    # Dimensions but not checked for consistency
    NDIM = ndims(mask)
    dim = size(mask)
    # Array which will hold the multiplication factors to be applied to Length scales
    lambda = ones(Float64, dim)
    dens[dens.<0] .= 0.
    dens = dens.^(1.0 / NDIM)

    dens[isnan.(dens)] .= 0
    # Total volume
    totalvolume = DIVAnd_integral(mask, pmn, lambda)
    #@show totalvolume
    denspower = DIVAnd_integral(mask, pmn, dens)
    #@show denspower
    lambda = (denspower / totalvolume) ./ dens
    lambda[.!mask] .= 1.0
    lambda[dens.==0.0] .= 1.0
    return lambda

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


