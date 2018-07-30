"""
Computes an estimate of the mean value of the diagonal of HK using GCV and the already solved analysis and it structure s

Only using real data locations


Kii = DIVAnd_GCVKiiobs(s);

"""
function DIVAnd_GCVKiiobs(s,nr=30;FIELD=())

    #the second, optional argument is the number of random vectors nr used for the estimate


    H = s.obsconstrain.H;
    R = s.obsconstrain.R;

	# if nr <0 use the data and analysis itself as random vector and KZ. Allows to use the function when s.P is not available
	#

	if nr<0
	nrealdata=sum(1 .- s.obsout);
    ndata=size(s.obsout)[1];

	#@show nrealdata,ndata,size(s.obsconstrain.yo)

	#@show size(s.yo),size(((s.H)*statevector_pack(s.sv,(FIELD,))))
	Kii=s.yo'*((s.H)*statevector_pack(s.sv,(FIELD,)))/(s.yo'*s.yo)

    if nrealdata==0
        Kii=0.0;
    else
        factorc=ndata/nrealdata;
        # Now take average of the nr different estimates,
        Kii=factorc*Kii;
    end
	return Kii

	end


    #if optimisation is to be used, make sure to use the same reference random points
    if VERSION >= v"0.7.0-beta.0"
        Random.seed!(nr)
    else
        srand(nr)
    end

    Z=randn(size(R)[1],nr);

    if VERSION >= v"0.7.0-beta.0"
        Random.seed!()
    else
        srand()
    end


    P = s.P;
    WW=P * (H'* (R \ Z));
    ZtHKZ=  Z'*(H*WW);
    ZtZ  =  Z'*Z;

    # correction for points out of the domain:
    nrealdata=sum(1 .- s.obsout);
    ndata=size(s.obsout)[1];
    if nrealdata==0
        Kii=0.0;
    else
        factorc=ndata/nrealdata;
        # Now take average of the nr different estimates,
        Kii=factorc*mean(diag(ZtHKZ)./diag(ZtZ));
    end
    return Kii

end

# Copyright (C) 2008-2017 Alexander Barth <barth.alexander@gmail.com>
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
