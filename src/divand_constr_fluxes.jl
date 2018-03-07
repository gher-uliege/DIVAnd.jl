# Creates integral constraints for each latitude so that a barotropic correction step leads 
# to an additional flux prescribed. 
#
# c = divand_constr_fluxes(s,topographyforfluxes,fluxes,epsfluxes,pmnin)
#
# Input:
#   s: structure
#   topographyforfluxes: 2D array with the bottom topography used for the flux calculations 
#               DO NOT USE NaN in it. 
#   fluxes: array of fluxes. The barotropic correction on elevation should be such that 
#                         Sum over longitude at each latidute of Sum h \delta(eta)/\delta x   \delta x = fluxes
#   epsfluxes: error variance on constraint. Scaling to be verified
#   pmnin: metrics from the calling routine
#
# 
# Output:
#   c: structure to be used by divand_addc with the following fields: R (a
#     covariance matrix), H (extraction operator) and yo (specified value for
#     the constrain).

function divand_constr_fluxes(s,topographyforfluxes,fluxes,epsfluxes,pmnin)

    

    mask = s.mask;

    n  = s.n;
    iscyclic = s.iscyclic;

    sz = size(mask);

# JMB directly put the right number of constraints for the moment, one for each latitude. sz[2] in the present case
# hardwired loop over all latitudes this is not a problem as long as only one direction is forced to have zero fluxes.
# Just rotate the domain before calculations.     

    jmjmax=sz[2]
	
#	A = spzeros(s.sv.n,s.sv.n)
    A = spzeros(jmjmax,s.sv.n)
	 l = size(A,1);
	 yo = zeros(l)
	 R = Diagonal(mean(topographyforfluxes)*epsfluxes.*ones(l))

    for j=1:jmjmax

    for i=1:1
        S = sparse_stagger(sz,i,iscyclic[i]);
        m = (S * mask[:]) .== 1;

        d = topographyforfluxes
# JMB: Add here integrals by using pack of an array with dx at a given latitude
# Take same shape as velocity array
        forintegral=zeros(d)
# Use metrics
        forintegral[:,j]=1.0./pmnin[1][:,j]
#  Pack forintegral
#
#        A = A + sparse_diag(d[mask]) * sparse_pack(mask) * S' * sparse_pack(m)' * s.Dx[i];
         packedline=statevector_pack(s.sv,(forintegral,))
		
		 jmw=packedline'*sparse_diag(d[mask]) * sparse_pack(mask) * S' * sparse_pack(m)' * s.Dx[i] 
		 #@show jmw
		 #@show size(jmw),size(A),size(A[j,:]),size(squeeze(jmw,1))
		 
		 
         A[j,:] = A[j,:] + squeeze(jmw,1);

# Simple test with fake value		 
         #yo[j]=j*100*100
    end
   yo=fluxes
   

	end
# end loop o j
	
    H = A;
    
    
    #R = speye(size(H,1));

    return divand_constrain(yo,R,H)

end

# Copyright (C) 2014, 2017 Alexander Barth <a.barth@ulg.ac.be>
#                     2018 Jean-Marie Beckers <JM.Beckers@uliege.be>
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
