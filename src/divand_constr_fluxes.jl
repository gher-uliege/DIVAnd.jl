# Creates integral constraints for each latitude so that a barotropic correction step leads
# to an additional flux prescribed.
#
# c = divand_constr_fluxes(s,topographyforfluxes,fluxes,epsfluxes,pmnin)
#
# Input:
#   s: structure
#   topographyforfluxes: tuple of two 2D arrays with the bottom topography used for the flux calculations
#               DO NOT USE NaN in it. If an array is replaced by a scalar zero, the constraint is not used.
#               for fluxes calculated with geostrophy apply g/f to h
#   fluxes: tuple of two arrays of fluxes. The barotropic correction on elevation should be such that
#                         Sum over longitude at each latidute of Sum h \delta(eta)/\delta x   \delta x = - fluxes[1]
#                         Sum over latitude  at each longitude of Sum h \delta(eta)/\delta y   \delta y = -fluxes[2]
#             WARNING: This has been coded to directly use geostrophy.jl output and flux directions
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

# JMB directly put the right number of constraints for the moment, one for each latitude.

    jmjmax=[0,0]
	if topographyforfluxes[1]!=0
	jmjmax[1]=sz[2]
	end

	if topographyforfluxes[2]!=0
	jmjmax[2]=sz[1]
	end

	if sum(jmjmax)==0
		warning("no constraint in _fluxes")
		return 0
	end


#	A = spzeros(s.sv.n,s.sv.n)
    A = spzeros(sum(jmjmax),s.sv.n)
	 l = size(A,1);
	 yo = zeros(l)
	 R = Diagonal((mean(topographyforfluxes[1])+mean(topographyforfluxes[2]))*epsfluxes.*ones(l))


     joffset=0
    for i=1:2
        S = sparse_stagger(sz,i,iscyclic[i]);
        m = (S * mask[:]) .== 1;

        d = topographyforfluxes[i]


	   for j=1:jmjmax[i]
# JMB: Add here integrals by using pack of an array with dx at a given latitude
# Take same shape as topo array
        forintegral=zeros(d)
# Use metrics
        if i==1
           forintegral[:,j]=1.0./pmnin[1][:,j]
		else
		   forintegral[j,:]=1.0./pmnin[2][j,:]
		end
#  Pack forintegral
#
#        A = A + sparse_diag(d[mask]) * sparse_pack(mask) * S' * sparse_pack(m)' * s.Dx[i];
         packedline=statevector_pack(s.sv,(forintegral,))

		 jmw=packedline'*sparse_diag(d[mask]) * sparse_pack(mask) * S' * sparse_pack(m)' * s.Dx[i]
		 #@show sparse_diag(d[mask]) * sparse_pack(mask) * S' * sparse_pack(m)' * s.Dx[i]
		 #@show mean(sparse_diag(d[mask]))
		 #@show size(jmw),size(A),size(A[j,:]),size(squeeze(jmw,1))

		 #@show squeeze(jmw,1),fluxes[i][j]

         A[j+joffset,:] = A[j+joffset,:] + squeeze(jmw,1);
		 # test is kept in case flux signs are changed to x,y instead normal direction
		 if i==1
         yo[j+joffset]=-fluxes[i][j]
		 else
		 yo[j+joffset]=-fluxes[i][j]
		 end

        end


		joffset=joffset+jmjmax[i]




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
