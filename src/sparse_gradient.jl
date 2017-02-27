"""
Sparse operator for a gradient.

Dx1,Dx2,...,Dxn = sparse_gradient(mask,pmn)

Form the gradient using finite differences in all n-dimensions

Input:
  mask: binary mask delimiting the domain. 1 is inside and 0 outside.
        For oceanographic application, this is the land-sea mask.

  pmn: scale factor of the grid.

Output:
  Dx1,Dx2,...,Dxn: operators represeting a gradient along
    different dimensions
"""
#JM added arguments
function sparse_gradient(operatortype,mask,pmn,Labs,alphabc,iscyclic = falses(ndims(mask)))

    H = oper_pack(operatortype,mask)

    sz = size(mask)
    n = ndims(mask)

    out = []

    for i=1:n
        # staggering operator
        S = oper_stagger(operatortype,sz,i,iscyclic[i])

        # mask for staggered variable
        m = (S * mask[:]) .== 1
#JM 
     if alphabc>0
#        d = m .* (S * pmn[i][:])

# Deep copy needed here otherwise changes up to the divandrun calling routine on pmn
# But on the other hand it seems to have other beneficial effects; so there must be another use of pmn somewhere ? In the laplacian near the border ?
         #wjmb=deepcopy(pmn[i])
		 wjmb=pmn[i]
# For the moment, hardcoded for 1D and 2D
	
#	    @show alphabc
        if n==1
		  if ~iscyclic[1]
            wjmb[1]=1.0/(alphabc.*Labs[1][1])
            wjmb[end]=1.0/(alphabc.*Labs[1][end])
		  end
#		  @show wjmb[1], pmn[1][1]
        end

        if n==2
		  if i==1
		  if ~iscyclic[1]
		    wjmb[1,:]=1.0./(alphabc.*Labs[1][1,:])
            wjmb[end,:]=1.0./(alphabc.*Labs[1][end,:])
		  end
		  end
		  if i==2
		  if ~iscyclic[2]
            wjmb[:,1]=1.0./(alphabc.*Labs[2][:,1])
            wjmb[:,end]=1.0./(alphabc.*Labs[2][:,end])
		  end
		  end
        end 
		d = m .* (S * wjmb[:])
	 else
	    d = m .* (S * pmn[i][:])
     end

       
#/JM
#

        push!(out,oper_pack(operatortype,m) * oper_diag(operatortype,d) * oper_diff(operatortype,sz,i,iscyclic[i]) * H')
    end

    (out...)

end


# Copyright (C) 2009,2016 Alexander Barth <a.barth@ulg.ac.be>
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; If not, see <http://www.gnu.org/licenses/>.
