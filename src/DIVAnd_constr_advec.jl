# Create the advection constrain.
#
# c = DIVAnd_constr_advec(s,velocity)
#
# Create the advection constrain using the specified velocity.
#
# Input:
#   s: structure created by DIVAnd_background
#   velocity: tuple of velocity vectors
#
# Output:
#   c: structure to be used by DIVAnd_addc with the following fields: R (a
#     covariance matrix), H (extraction operator) and yo (specified value for
#     the constrain).

function DIVAnd_constr_advec(s,velocity)

    # check for NaNs
    nancount = sum([sum(isnan.(v)) for v in velocity])
    if nancount > 0
        error("$(nancount) velocity values are equal to NaN");
    end

    mask = s.mask;

    n  = s.n;
    iscyclic = s.iscyclic;

    sz = size(mask);

    A = spzeros(s.sv.n,s.sv.n)

    for i=1:n
        S = sparse_stagger(sz,i,iscyclic[i]);
        m = (S * mask[:]) .== 1;

        d = velocity[i]

        A = A + sparse_diag(d[mask]) * sparse_pack(mask) * S' * sparse_pack(m)' * s.Dx[i];
    end

    l = size(A,1);

    H = A;
    yo = zeros(l)
    R = Diagonal(ones(l))
    #R = speye(size(H,1));

    return DIVAnd_constrain(yo,R,H)

end

# Copyright (C) 2014, 2017 Alexander Barth <a.barth@ulg.ac.be>
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
