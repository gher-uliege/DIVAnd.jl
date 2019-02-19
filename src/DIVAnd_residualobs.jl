"""
    dataresidual = DIVAnd_residualobs(s,fi);

Computes the residual yo - H xa  only at real data points using the analysis.
on the grid `fi` and the solution structure `s`.
"""
function DIVAnd_residualobs(s,fi)
    residual = zeros(length(s.obsconstrain.yo))
    residual .= s.obsconstrain.yo-(s.obsconstrain.H)*statevector_pack(s.sv,(fi,))
    residual[s.obsout] .= NaN
    return residual
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
