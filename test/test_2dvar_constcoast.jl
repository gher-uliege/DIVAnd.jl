# Testing divand in 2 dimensions with advection.
using divand
using Base.Test

# grid of background field
mask,(pm,pn),(xi,yi) = divand_squaredom(2,linspace(-1,1,30))

# island at these location
mi0 = 12
mi1 = 13

mj0 = 12
mj1 = 13

mask[mi0:mi1,mj0:mj1] = false

x = [0.]
y = [0.]
f = [1.]

epsilon2 = 1/200
len = 0.3

fi,s = divandrun(mask,(pm,pn),(xi,yi),(x,y),f,len,epsilon2);


#=

The gradient in drection x at the point > should be close to zero if one 
of the grid cell denoted by * is a land point.

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
  |                 |        .        >                 |                 |
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

function boundaryu(mask)

    bu = falses(size(mask,1)-1,size(mask,2))
    
    for j = 1:size(bu,2)
        for i = 1:size(bu,1)
            # check if i+½,j is a see point
            
            if mask[i,j] && mask[i+1,j]
            
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
        
    end
    return bu
end

boundaryv(mask) = boundaryu(mask')'


function normalvel(mask,fip)
    fi = zeros(size(mask))
    fi[mask] = fip
    bu = boundaryu(mask)
    bv = boundaryv(mask)
    Dxfi = fi[2:end,:] - fi[1:end-1,:]
    Dyfi = fi[:,2:end] - fi[:,1:end-1]

    return [Dxfi[bu]; Dyfi[bv]]
end


"""
    H = normalvelsp(mask)

Return a sparse matrix `H` computing the
gradient of the values along the coastline defined by `mask`.
Land points correspond to the value `false` and sea points to the value `true`
in `mask`.
"""
function normalvelsp(mask)
    sz = size(mask)
    Hunpack = divand.sparse_pack(mask)'
    bu = boundaryu(mask)
    bv = boundaryv(mask)

    diffx = sparse_diff(sz,1)
    diffy = sparse_diff(sz,2)
    
    H = [divand.sparse_pack(bu) * diffx * Hunpack;
         divand.sparse_pack(bv) * diffy * Hunpack]
    
    return H
end

function divand_constcoast(mask,eps2)
    H = normalvelsp(mask)
    m = size(H,1)
    yo = zeros(m)
    R = Diagonal(fill(eps2,(m,)))
    
    return divand.divand_constrain(yo,R,H)
end

eps2 = 1e-7

c = divand_constcoast(mask,eps2)

@test normalvel(mask,fi[mask]) ≈ normalvelsp(mask) * fi[mask]


fi2,s = divandrun(mask,(pm,pn),(xi,yi),(x,y),f,len,epsilon2;
                  constraints = [c]);


# extract boundary values around the island and check
# if the standard deviation is small

@test std(fi2[mi0-1:mi1+1,mj0-1:mj1+1][mask[mi0-1:mi1+1,mj0-1:mj1+1]]) < 1e-5


# more complex example

srand(1234)
mask0,(pm,pn),(xi,yi) = divand_squaredom(2,linspace(0,1,100))
mask = divand.random(mask0,(pm,pn),0.1,1)[:,:,1] .> 0.5
x = rand(100)
y = rand(size(x))
f = sin.(2*π*x) .* sin.(2*π*y)

len = 0.1

fi3,s = divandrun(mask,(pm,pn),(xi,yi),(x,y),f,len,epsilon2;
                  constraints = [divand_constcoast(mask,eps2)]);


nothing



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
