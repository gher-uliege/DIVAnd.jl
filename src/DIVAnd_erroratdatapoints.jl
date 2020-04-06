"""
    errorvariance = DIVAnd_erroratdatapoints(s;restrictedlist=[])

Computes the error at the real data locations using the analysis structure s

If a restricedlist is provided erros are only calculated at the indexes where restricedlist==true

"""
function DIVAnd_erroratdatapoints(s; restrictedlist = [])
    if restrictedlist == []
        return diagMtCM(s.P, s.obsconstrain.H')
    else
        ei = zeros(Float64, (size(restrictedlist)[1], 1))
        errorat = zeros(Float64, size(restrictedlist)[1])
        for i in Iterators.filter(x -> x[2] == true, enumerate(restrictedlist))
            ei[:, 1] .= 0
            ei[i[1], 1] = 1
            errorat[i[1]] = diagMtCM(s.P, s.obsconstrain.H' * ei)[1]
        end
        return errorat
    end
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

# LocalWords:  fi DIVAnd pmn len diag CovarParam vel ceil moddim fracdim
