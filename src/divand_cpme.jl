"""
Compute a variational analysis of arbitrarily located observations to calculate the clever poor man's error

cpme = divand_cpme(mask,pmn,xi,x,f,len,lambda,...);

Perform an n-dimensional variational analysis of the observations `f` located at
the coordinates `x`. The array `cpme` represent the error field at the grid
defined by the coordinates `xi` and the scales factors `pmn`.

# Input:
* `mask`: binary mask delimiting the domain. true is inside and false outside. For oceanographic application, this is the land-sea mask.

* `pmn`: scale factor of the grid. pmn is a tuple with n elements. Every
       element represents the scale factor of the corresponding dimension. Its
       inverse is the local resolution of the grid in a particular dimension.

*  `xi`: tuple with n elements. Every element represents a coordinate
  of the final grid on which the observations are interpolated

* `x`: tuple with n elements. Every element represents a coordinate of
  the observations

* `f`: value of the observations *minus* the background estimate (m-by-1 array).
    (see note)

* `len`: correlation length

* `lambda`: signal-to-noise ratio of observations (if lambda is a scalar).
    The larger this value is, the closer is the field `fi` to the
    observation. If lambda is a scalar, then R is 1/lambda I, where R is the observation error covariance matrix). If lambda is a vector, then R is diag(lambda) or if lambda is a matrix (a matrix-like project), then R is equal to lambda.


# Optional input arguments specified as keyword arguments also as for divand
"""


function divand_cpme(mask,pmn,xi,x,f,len,lambda; otherargs...)

# check inputs

if !any(mask[:])
  error("no sea points in mask");
end

errorscale=1;

cpme,s =  divandrun(mask,pmn,xi,x,ones(size(f)),len./1.70766,lambda; otherargs...);
cpme=errorscale.*(-cpme.+1);

return cpme

end

# Copyright (C) 2008-2017 Alexander Barth <barth.alexander@gmail.com>
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

# LocalWords:  fi divand pmn len diag CovarParam vel ceil moddim fracdim


