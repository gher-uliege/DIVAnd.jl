"""
   s = DIVAnd_addc(s,c)

Add a constraint `c` to the cost function defined by `s`.
The structure `s` is typically created by DIVAnd_background and the contrain `c`
has the following fields: `R` (a covariance matrix), `H` (extraction operator) and
`yo` (specified value for the constrain).
The added contrain Jc(x) is quadratic and has the following structure.

Jc(x) = (H x - yo)ᵀ R⁻¹ (H x - yo)

"""
function DIVAnd_addc(s,constrain)
    if isempty(s.H)
        s.H = constrain.H
        s.R = constrain.R
        s.yo = constrain.yo
    else
        s.H = vcat(s.H,constrain.H)
        s.R = blkdiag(s.R,constrain.R)
        s.yo = vcat(s.yo,constrain.yo)
    end
    return s
end

# Copyright (C) 2014,2019 Alexander Barth <a.barth@ulg.ac.be>
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
