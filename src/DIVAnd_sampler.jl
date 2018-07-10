"""


samplesteps = DIVAnd_sampler(pmn,len);

# Defines steps for sub-sampling in the discrete grid which would still allow to resolve the provided lengthscales

# Input:

* `pmn`: scale factor of the grid. pmn is a tuple with n elements. Every
       element represents the scale factor of the corresponding dimension. Its
       inverse is the local resolution of the grid in a particular dimension.

* `len`: correlation length



# Output:

* `samplesteps`: vector of integers with steps in subsampling [1 2 4 1] means every grid point in x direction, every fifth in y etc

"""
function DIVAnd_sampler(pmn,len)

    # TO DO: in a single sweep compute both minimum and maximum with function extrema
    # Usefull for DIVAndgo and DIVAndjog


    n = ndims(pmn[1])
    samplesteps=ones(Int,n);
    Labspmnmin=zeros(n)

    for i=1:n
        if isa(len,Number)
            Labspmnmin[i] = len*minimum(pmn[i]);
        elseif isa(len,Tuple)

            if isa(len[1],Number)
                Labspmnmin[i] = len[i]*minimum(pmn[i]);

            else
                Labspmnmin[i] = minimum(len[i].*pmn[i])

            end

        end

        nsamp=Int(floor(Labspmnmin[i]/3));
        if nsamp>1
            samplesteps[i]=nsamp
        end
    end



    return samplesteps



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
