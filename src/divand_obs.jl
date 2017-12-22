"""
Include the constrain from the observations.

s = divand_obs(s,xi,x,R,I)

Set observations of variational problem.
It is assumed that the each coordinate depends only on one
index. If this is not the case, then matrix I must be provided.

Input:
  s: structure created by divand_background
  xi: coordinates of observations*
  x: coordinates of grid*
  R: obs. error covariance matrix (normalized)
  I (optional): fractional indexes of location of observation
    within the grid

Output:
  s: structure to be used by divand_factorize

Note:
  *these parameters can either be specified as a cell
  array of all dimenions:
  xi = {Xi,Yi,Zi}
  or as n+1 dimensional array
"""
function divand_obs(s,xi,x,yo,R; I = [])


    #xi = cat_cell_array(xi);
    #x = cat_cell_array(x);

    mask = s.mask;
    iscyclic = s.iscyclic;
    moddim = s.moddim;


    if isempty(I)
        I = localize_separable_grid(x,mask,xi);
    end

    H,out,outbbox = sparse_interp(mask,I,iscyclic);

    nout = sum(out);
    if nout != 0
        noutbbox = sum(outbbox);
        #warn("Observations out of bounding box: $(noutbbox) and touching land $(nout-noutbbox)")
    end

    # NaN points
    nanobs = isnan.(yo)
    nnanobs = sum(nanobs)
    if nnanobs != 0
        out = out .| nanobs
        yo = deepcopy(yo)
        yo[nanobs] = 0.
        warn("Observations equal to NaN: $(nnanobs)")
    end

    H = H * sparse_pack(mask)';


    s.obsout = out;
    constrain = divand_constrain(yo,R,H)
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
