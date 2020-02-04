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
function DIVAnd_adaptedeps2(s, fi)

    residual = DIVAnd_residualobs(s, fi)
    diagR = diag(s.obsconstrain.R)::Vector{Float64}

    return DIVAnd_adaptedeps2(s.yo, residual, diagR, s.obsout)
end


"""
    DIVAnd_adaptedeps2(yo, residual, diagR, ignoreobs)

Using Deroziers adaptive approach provides a multiplicative factor for the current epsilon2 value so that factor*epsilon2 is a better
estimate of the R matrix.

`yo` the observations (minus the background),
`residual` the obserations minus the analysis,
`diagR`, the diagonal of the rel. obs. error covariance matrix and
`ignoreobs` is true if an observation is out of the grid or should be ignored for other reasons.

For unscaled R and assuming that the background is zero, Deroziers showed that:

mean((yᵒ - Hxᵃ) ⋅ yᵒ) =  ϵ²
mean(yᵒ ⋅ yᵒ) = σ² +  ϵ²

mean(yᵒ ⋅ yᵒ) / mean((yᵒ - Hxᵃ) ⋅ yᵒ) = σ²/ϵ² + 1
λ = σ²/ϵ² = 1 - mean(yᵒ ⋅ yᵒ) / mean((yᵒ - Hxᵃ) ⋅ yᵒ)

ϵ² / σ² = 1 / λ
"""
function DIVAnd_adaptedeps2(yo, residual, diagR, ignoreobs)
    d0d = zero(eltype(yo))     # yo ⋅ yo
    d0dmd1d = zero(eltype(yo)) # (yo - Hxa) ⋅ yo
    eps2 = zero(eltype(yo))    # mean of diagonal of scaled R
    inv_eps2 = zero(eltype(yo))    # mean of inv of the diagonal of scaled R
    nrealdata = 0

    for i = 1:length(yo)
        if !ignoreobs[i]
            d0d += yo[i]^2 / diagR[i]
            d0dmd1d += yo[i] * residual[i] / diagR[i]
            inv_eps2 += 1/diagR[i]
            nrealdata += 1
        end
    end

    inv_eps2 /= nrealdata

    ll1 = d0d / (d0dmd1d) - 1
    eps1 = 1 / ll1
    factor = eps1 * inv_eps2


    if (factor == 0) || !isfinite(factor)
        error("scalefactore has the value $factor")
    end

    return factor
end


DIVAnd_adaptedeps2(yo, residual, diagR::Number, ignoreobs) =
    DIVAnd_adaptedeps2(yo, residual, fill(diagR,size(yo)), ignoreobs)

DIVAnd_adaptedeps2(yo, residual, R::Matrix, ignoreobs) =
    DIVAnd_adaptedeps2(yo, residual, diag(R), ignoreobs)

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
