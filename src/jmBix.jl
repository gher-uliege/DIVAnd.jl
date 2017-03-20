"""

"""

function jmBix(s,x;btrunc=[])

    iBx=s.iB*x
    # @show btrunc
    if btrunc==[]
        return iBx
    end


    WE = s.WE;
    coeff = s.coeff;
    n = s.n;
    alpha=s.alpha
    # incomplete Bi calculated before now complemented on the fly

    D=s.D
    for j=btrunc+1:length(alpha)
        # exponent of laplacian
        k = Int(floor((j-2)/2))
        # @show size(iBx)
        # @show size(x)
        # @show size(s.iB)
        #iBx_ = spzeros(iBx);
        iBx_=0.*iBx;

        if mod(j,2) == 0
            # constrain of derivative with uneven order (j-1)
            # (gradient, gradient*laplacian,...)
            # normalized by surface
            Dk=D^k
            Dkx=Dk*x
            for i=1:n
                #                Dx = s.WEss[i] * (s.Dx[i] * Dk);
                #                iBx_ = iBx_ + Dx'*(Dx*x);
                Dx = s.WEss[i] * s.Dx[i];
                # maybe gain if Dk=Dk'?
                iBx_ = iBx_ + Dk'*(Dx'*(Dx*Dkx));
            end

        else
            # constrain of derivative with even order (j-1)
            # (laplacian, biharmonic,...)

            # normalize by surface of each cell
            # such that inner produces (i.e. WE'*WE)
            # become integrals
            # WD: units length^(n/2)

            WD = WE * D^(k+1);
            iBx_ = WD'*(WD*x);
        end

        iBx_ = iBx_/coeff;



        iBx = iBx + alpha[j] * iBx_
    end
    return iBx
end

# LocalWords:  iB divand

# Copyright (C) 2014 Alexander Barth <a.barth@ulg.ac.be>
#                         Jean-Marie Beckers   <JM.Beckers@ulg.ac.be>
#
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
