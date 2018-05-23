# Testing divand in 2 dimensions with advection.
using divand
using Base.Test

# grid of background field
mask,(pm,pn),(xi,yi) = divand_squaredom(2,linspace(-1,1,30))
mask[5:10,5:10] = false


#mask,(pm,pn),(xi,yi) = divand_squaredom(2,linspace(-1,1,3))
#mask[2,2:3] = false
#mask[3,2] = false

x = [.4]
y = [.4]
f = [1.]

a = 5;
u = a*yi;
v = -a*xi;
epsilon2 = 1/200
len = 0.2

fi,s = divandrun(mask,(pm,pn),(xi,yi),(x,y),f,len,epsilon2;velocity = (u,v),alphabc=0);

Dxfi = fi[2:end,:] - fi[1:end-1,:]
Dyfi = fi[:,2:end] - fi[:,1:end-1]

#boundary_u, boundary_v, boundary_psi = stagger_mask(mask,xor)
#Dxfi[

function boundaryu(mask)

    bu = falses(size(mask,1)-1,size(mask,2))
    
    for j = 1:size(bu,2)
        for i = 1:size(bu,1)
            
            if j > 1
                if !mask[i,j-1] || !mask[i+1,j-1]
                    bu[i,j] = true
                end
            end
            
            if j < size(bu,2)
                if !mask[i,j+1] || !mask[i+1,j+1]
                    bu[i,j] = true
                end            
            end
        end
        
    end
    return bu
end

bu = boundaryu(mask)
bv = boundaryu(mask')'

#=

  |                 |                 |                 |                 |
  +-----------------+-----------------+-----------------+-----------------+
  |                 |                 |                 |                 |
  |                 |                 |                 |                 |
  |                 |                 |                 |                 |
  |                 |        *        |        *        |                 |
  |                 |     (i,j+1)     |    (i+1,j+1)    |                 |
  |                 |                 |                 |                 |
  |                 |                 |                 |                 |
  +-----------------+--------^--------+-----------------+-----------------+
  |                 |      (i,j)      |                 |                 |
  |                 |                 |                 |                 |
  |                 |                 |                 |                 |
  |                 |        *        >                 |                 |
  |                 |      (i,j)    (i,j)               |                 |
  |                 |                 |                 |                 |
  |                 |                 |                 |                 |
  +-----------------+-----------------+-----------------+-----------------+
  |                 |                 |                 |                 |
  |                 |                 |                 |                 |
  |                 |                 |                 |                 |
  |                 |        *        |        *        |                 |
  |                 |     (i,j-1)     |    (i+1,j-1)    |                 |
  |                 |                 |                 |                 |
  |                 |                 |                 |                 |
  +-----------------+-----------------+-----------------+-----------------+






=#
function normalvel(mask,fip)
    fi = zeros(size(mask))
    fi[mask] = fip
    bu = boundaryu(mask)
    bv = boundaryu(mask')'
    Dxfi = fi[2:end,:] - fi[1:end-1,:]
    Dyfi = fi[:,2:end] - fi[:,1:end-1]

    return [Dxfi[bu]; Dyfi[bv]]
end


function normalvelsp(mask,fip)
    sz = size(mask)
    fi = zeros(size(mask))
    fi[mask] = fip

    Hunpack = divand.sparse_pack(mask)'
    fi = reshape(divand.sparse_pack(mask)' * fip,size(mask))
    bu = boundaryu(mask)
    bv = boundaryu(mask')'

    divand.sparse_pack(mask)
    diffx = sparse_diff(sz,1)
    diffy = sparse_diff(sz,2)
    
    Dxfi = fi[2:end,:] - fi[1:end-1,:]
    @show size(Dxfi)
    
    Dxfi2 = divand.sparse_pack(bu) * diffx * Hunpack * fip
    Dyfi2 = divand.sparse_pack(bv) * diffy * Hunpack * fip
    @show size(Dxfi2)
    Dyfi = fi[:,2:end] - fi[:,1:end-1]


    H = [divand.sparse_pack(bu) * diffx * Hunpack;
         divand.sparse_pack(bv) * diffy * Hunpack]
    
    #return [Dxfi2; Dyfi2]
    return H*fip
end

@test normalvel(mask,fi[mask]) ≈ normalvelsp(mask,fi[mask])


#nothing
#@test abs(fi[18,24] - 0.8993529043140029) < 1e-2



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
