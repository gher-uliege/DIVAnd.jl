# Form the different components of the background error covariance matrix.
#
# [iB_,iB] = divand_background_components(s,alpha)
#
# Compute the components of the background error covariance matrix iB_ and
# their sum based on alpha (the a-dimensional coefficients for norm, gradient,
# laplacian,...).

function divand_background_components(s,alpha)

WE = s.WE;
D = s.D;
coeff = s.coeff;
n = s.n;

iB_ = Array{SparseMatrixCSC{Float64,Int64}}(length(alpha));

# constrain of total norm

iB_[1] =  (1/coeff) * (WE'*WE);


# loop over all derivatives

for j=2:length(alpha)
    # exponent of laplacian
    k = Int(floor((j-2)/2))

    if mod(j,2) == 0
        # constrain of derivative with uneven order (j-1)
        # (gradient, gradient*laplacian,...)
        # normalized by surface

        iB_[j] = spzeros(size(D,1),size(D,1));

        for i=1:n
            Dx = s.WEss[i] * s.Dx[i] * D^k;
            #Dx = s.WEs[i] * s.Dxs[i] * D^k;
            iB_[j] = iB_[j] + Dx'*Dx;
        end
    else
        # constrain of derivative with even order (j-1)
        # (laplacian, biharmonic,...)

        # normalize by surface of each cell
        # such that inner produces (i.e. WE'*WE)
        # become integrals
        # WD: units length^(n/2)

        WD = WE * D^(k+1);
        iB_[j] = WD'*WD;
    end

    iB_[j] = iB_[j]/coeff;
end


# sum all terms of iB
# iB is adimentional
iB = alpha[1] * iB_[1]
for j=2:length(alpha)
  iB = iB + alpha[j] * iB_[j]
end


return iB_,iB
end

# LocalWords:  iB divand

# Copyright (C) 2014 Alexander Barth <a.barth@ulg.ac.be>
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
