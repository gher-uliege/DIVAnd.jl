"""
    s = DIVAnd_obs(s,xi,x,R,I)

Include the constrain from the observations.
It is assumed that the each coordinate depends only on one
index. If this is not the case, then matrix I must be provided.

Input:
  s: structure created by DIVAnd_background
  xi: coordinates of observations (tuple of vectors)
  x: coordinates of grid (tuple of arrays)
  R: obs. error covariance matrix (normalized)
  I (optional): fractional indexes of location of observation
    within the grid

Output:
  s: structure to be used by DIVAnd_factorize

Note make sure not to mix Float32 and Float64 for DIVAnd_constrain.
"""
function DIVAnd_obs(s,xi,x,yo::Vector{T},R,I = zeros(T,0,0)) where T
    mask = s.mask
    iscyclic = s.iscyclic
    moddim = s.moddim

    if isempty(I)
        I = localize_separable_grid(x,mask,xi) :: Matrix{T}
    end

    H,out,outbbox = sparse_interp(mask,I,iscyclic)

    nout = sum(out)
    if nout != 0
        noutbbox = sum(outbbox)
        #@warn "Observations out of bounding box: $(noutbbox) and touching land $(nout-noutbbox)"
    end

    # NaN points
    nanobs = isnan.(yo)
    nnanobs = sum(nanobs)
    if nnanobs != 0
        out = out .| nanobs
        yo = deepcopy(yo)
        yo[nanobs] .= 0.
        @warn "Observations equal to NaN: $(nnanobs)"
    end

    H = H * sparse_pack(mask)'

    s.obsout = out

    if isa(R,Diagonal)
        diagR = Float64.(diag(R))
        diagR[out] .= Inf
        R = Diagonal(diagR)
    else
        error("all observation must be inside the domain for non-diagonal error observation covariance matrix")
    end

    constrain = DIVAnd_constrain(yo,R,H)

    s.obsconstrain = constrain

    return constrain
end

# Copyright (C) 2014,2017 Alexander Barth <a.barth@ulg.ac.be>
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
