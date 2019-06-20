"""

    factor = DIVAnd_adaptedeps2(s,fi);

# Input:
* `s`: structure returned by `DIVAndrun`
* `fi`: analysis returned by `DIVAndrun`

# Output:
* `factor` : multiplicative factor to apply to epsilon2


Using Deroziers adaptive approach provides a multiplicative factor for the current epsilon2 value so that factor*epsilon2 is a better
estimate of the R matrix. If you cannot use `DIVAndrun` but use `DIVAndgo`, the latter provides automatically this pamater as result.


"""
function DIVAnd_adaptedeps2(s,fi);

    residual = DIVAnd_residualobs(s,fi);

    d0d = zero(eltype(fi))     # yo ⋅ yo
    d0dmd1d = zero(eltype(fi)) # (yo - Hxa) ⋅ yo
    diagR = diag(s.obsconstrain.R) :: Vector{Float64}
    eps2 = zero(eltype(fi))    # mean of diagonal of R
    nrealdata = 0

    for i = 1:length(s.yo)
        if !s.obsout[i]
            d0d +=  s.yo[i]^2
            d0dmd1d +=  s.yo[i] * residual[i]
            eps2 += diagR[i]
            nrealdata += 1
        end
    end
    eps2 /= nrealdata

    ll1 = d0d/(d0dmd1d)-1;
    eps1 = 1/ll1;
    factor = eps1/eps2;
    return factor
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
