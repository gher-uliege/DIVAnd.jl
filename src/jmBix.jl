function jmBix(s,x::Array{Float64,1};btrunc=[])
    iBx=s.iB*x   ::Array{Float64,1}

    if btrunc==[]
        return iBx
    end

    #WE = s.WE;
    coeff = s.coeff;
    n = s.n;
    alpha=s.alpha
    # incomplete Bi calculated before now complemented on the fly

    #D=s.D
    for j=btrunc+1:length(alpha)
        # exponent of laplacian
        k = Int(floor((j-2)/2))
        # @show size(iBx)
        # @show size(x)
        # @show size(s.iB)
        #iBx_ = spzeros(iBx);
        iBx_=0 * iBx ::Array{Float64,1};
        # Certainly a gain to make; not recompute D^(k+1) but Dk*D if k+1 is one larger than already calculated value if it exists
        # But only for n larger than 5 probably, so not an urgent thing

        if mod(j,2) == 0
            # constrain of derivative with uneven order (j-1)
            # (gradient, gradient*laplacian,...)
            # normalized by surface

			Dk=s.D^k
			Dkx=Dk*x ::Array{Float64,1}

            # see dirty hack in DIVAnd_background_components

            #=
            for i=1:n
			    if s.Ld[i]>0
				    iBx_ = iBx_ + (s.Dx[i]'*(s.WEss[i] *(s.WEss[i] *(s.Dx[i]*Dkx))));
				end
			end
            =#
			iBx_ = iBx_ + s.WEss[1]*Dkx

			if k>0
				iBx_=Dk'*iBx_
			end
        else
            # constrain of derivative with even order (j-1)
            # (laplacian, biharmonic,...)

            # normalize by surface of each cell
            # such that inner produces (i.e. WE'*WE)
            # become integrals
            # WD: units length^(n/2)
            #@show k
            #@time WD = s.WE * s.D^(k+1);
		    sDkp=s.D^(k+1)
			#iBx_ = WD'*(WD*x);

			iBx_ = sDkp'*(s.WE*(s.WE*(sDkp*x)))

			#Dk=s.WE * s.D^(k+1)
            #Dkx=Dk*x
            #iBx_ = Dk'*Dkx
        end
        asurc=alpha[j]/coeff
        #iBx_ = iBx_/coeff;
		#iBx = iBx + alpha[j] * iBx_

        iBx=BLAS.axpy!(asurc,iBx_,iBx)

    end
    return iBx
end

# LocalWords:  iB DIVAnd

# Copyright (C) 2014-2017 Alexander Barth	  	 <a.barth@ulg.ac.be>
#                         Jean-Marie Beckers 	 <JM.Beckers@ulg.ac.be>
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
