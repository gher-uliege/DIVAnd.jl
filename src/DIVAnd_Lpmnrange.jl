"""


Lpmnrange = DIVAnd_Lpmnrange(pmn,len);

# In each direction searches for the minimum and maximum value of the length scale times the metric in this diretion
# Si it basically looks at the worst and the best resolution found in the grid

# Input:

* `pmn`: scale factor of the grid. pmn is a tuple with n elements. Every
       element represents the scale factor of the corresponding dimension. Its
       inverse is the local resolution of the grid in a particular dimension.

* `len`: correlation length



# Output:

* `Lpmnrange`: Array of range tuples (minimum and maximum of L times metric)

"""
function DIVAnd_Lpmnrange(pmn::NTuple{N,Array{T,N}},len) where {N,T}
    Lpmnrange = Vector{NTuple{2,T}}(undef,N)

    for i=1:N
        if isa(len,Number)
            Lpmnrange[i] = extrema(len*pmn[i]);
        elseif isa(len,Tuple)

            if isa(len[1],Number)
                Lpmnrange[i] = extrema(len[i]*pmn[i]);
            else
                Lpmnrange[i] = extrema(len[i].*pmn[i])
            end

        end
    end

    return Lpmnrange
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
