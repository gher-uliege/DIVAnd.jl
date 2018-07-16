"""


    cpme = DIVAnd_cpme(mask,pmn,xi,x,f,len,epsilon2;...);



# Input: Same as for `DIVAndrun`
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

* `epsilon2`: error variance of the observations (normalized by the error variance of the background field). `epsilon2` can be a scalar (all observations have the same error variance and their errors are decorrelated), a vector (all observations can have a difference error variance and their errors are decorrelated) or a matrix (all observations can have a difference error variance and their errors can be correlated). If `epsilon2` is a scalar, it is thus the *inverse of the signal-to-noise ratio*.

* `keywords` : undocumented for the moment how to use iterative solver with coarser grid as preconditionner. see `DIVAndjog` for `csteps`, `lmask` and `alphapc`parameters

# Optional input arguments specified as keyword arguments also as for DIVAnd


# Output:

* `cpme`: the clever poor mans error

Perform an n-dimensional variational analysis of the observations `f` located at
the coordinates `x`. The array `cpme` represent the error field at the grid
defined by the coordinates `xi` and the scales factors `pmn`. If you cannot run `DIVAndrun` you can use `DIVAndgo` with error field calculation `:cpme`

"""
function DIVAnd_cpme(mask,pmn,xi,x,f,Labs,epsilon2; csteps=[0],lmask=[], alphapc=[], otherargs...)
#function DIVAnd_cpme(mask,pmn,xi,x,f,Labs,epsilon2)
#    csteps=[0]; lmask=[]; alphapc=[]; otherargs = Dict()

    errorscale=1;

    # The factor 1.70677 is the best one in 2D but should be slightly different for other dimensions
    # Could be a small improvement. Also used in DIVAnd_aexerr

    len = len_harmonize(Labs,mask)
    for i = 1:length(len)
        len[i] .= len[i] / 1.70766
    end

    if sum(csteps)>0
        cpme,s =  DIVAndjog(mask,pmn,xi,x,ones(size(f)),len,epsilon2,csteps,lmask; alphapc=alphapc, otherargs...);
    else
        cpme,s =  DIVAndrun(mask,pmn,xi,x,ones(size(f)),len,epsilon2; otherargs...);
    end
    cpme=errorscale .* max.(-cpme.+1,0)

    return cpme

end

# Copyright (C) 2008-2017 Alexander Barth <barth.alexander@gmail.com>
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

# LocalWords:  fi DIVAnd pmn len diag CovarParam vel ceil moddim fracdim
