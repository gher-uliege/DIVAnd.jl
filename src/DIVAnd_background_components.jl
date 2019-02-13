"""
    iB = DIVAnd_background_components(s,D,alpha; kwargs...)

Form the different components of the background error covariance matrix.
Compute the components of the background error covariance matrix `s.iB_` and
their sum based on alpha (the adimensional coefficients for norm, gradient,
laplacian,...).

If the optional arguments contains btrunc, the calculation of iB is limited
to the term up and including alpha[btrunc]
"""
function DIVAnd_background_components(s,D,alpha;
                                      coeff_derivative2::Vector{Float64} = zeros(ndims(mask)),
                                      kwargs...)

    WE = s.WE;
    coeff = s.coeff;
    n = s.n;

    kw = Dict(kwargs)

    # constrain of total norm

    iB_ =  (1/coeff) * (WE'*WE)

    if haskey(kw,:iB)
        kw[:iB][1] = iB_
    end

    # Truncate stored iB AFTER term btrunc, so alpha[btrunc] is the last one
    btrunc=length(alpha)
    if haskey(kw,:btrunc)
        btruncv=kw[:btrunc]
        if btruncv != Any[]
		    if btrunc>btruncv
				btrunc=btruncv
			end
        end
    end

    # sum all terms of iB
    # iB is adimentional
    iB = alpha[1] * iB_

    # loop over all derivatives

    for j=2:btrunc

        # exponent of laplacian
        k = Int(floor((j-2)/2))

        iB_ = spzeros(size(D,1),size(D,1));

        if mod(j,2) == 0
            # constrain of derivative with uneven order (j-1)
            # (gradient, gradient*laplacian,...)
            # normalized by surface

            for i=1:n
			# OPTIMIZATION: Do not calculate in directions where L is zero
			   if s.Ld[i]>0
                Dx = s.WEss[i] * s.Dx[i] * D^k;
                iB_ = iB_ + Dx'*Dx;
			   end
            end

        else
            # constrain of derivative with even order (j-1)
            # (laplacian, biharmonic,...)

            # normalize by surface of each cell
            # such that inner produces (i.e. WE'*WE)
            # become integrals
            # WD: units length^(n/2)

            WD = WE * D^(k+1);
            iB_ = WD'*WD;
        end

        iB_ = iB_/coeff;

        if haskey(kw,:iB)
            kw[:iB][j] = iB_
        end

        iB = iB + alpha[j] * iB_

    end

    # iB is adimensional

    #
# JMB Dirty hack 15/05: now s.WEss[1] contains the sum of all terms
    s.WEss[1]=s.Dx[1]'*(s.WEss[1] *(s.WEss[1] *(s.Dx[1])))
    for i=2:n
	  if s.Ld[i]>0
		s.WEss[1]=s.WEss[1]+s.Dx[i]'*(s.WEss[i] *(s.WEss[i] *(s.Dx[i])))
	  end
	end

    return iB
end

# LocalWords:  iB DIVAnd

# Copyright (C) 2014-2017 Alexander Barth	 <a.barth@ulg.ac.be>
#                         Jean-Marie Beckers <JM.Beckers@ulg.ac.be>
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
