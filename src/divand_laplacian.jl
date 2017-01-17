# Create the laplacian operator.
# 
# Lap = divand_laplacian(mask,pmn,nu,iscyclic)
#
# Form a Laplacian using finite differences
# assumes that gradient is zero at "coastline"
#
# Input:
#   mask: binary mask delimiting the domain. 1 is inside and 0 outside.
#         For oceanographic application, this is the land-sea mask.
#   pmn: scale factor of the grid. 
#   nu: diffusion coefficient of the Laplacian
#      field of the size mask or cell arrays of fields
#
# Output:
#   Lap: sparce matrix represeting a Laplaciant 
#
# 
function divand_laplacian(mask,pmn,nu,iscyclic)

# number of dimensions
n = ndims(mask)
sz = size(mask)

if !isa(nu,Tuple)
  # assume diffusion is the same along all dimensions
  nu = ([nu for i = 1:n]...)
end

# extraction operator of sea points
H = sparse_pack(mask)
sz = size(mask)

DD = sparse(Array{Int64}([]),Array{Int64}([]),Array{Float64}([]),prod(sz),prod(sz))

for i=1:n
  # operator for staggering in dimension i
  S = sparse_stagger(sz,i,iscyclic[i])

  # d = 1 for interfaces surounded by water, zero otherwise
  d = (S * mask[:]) .== 1
  
  # metric
  for j = 1:n
    tmp = S * pmn[j][:]
    
    if j==i
      d = d .* tmp
    else
      d = d ./ tmp
    end  
  end
  
  # nu[i] "diffusion coefficient"  along dimension i  

  d = d .* (S * nu[i][:])

  # Flux operators D
  # zero at "coastline"
  
  D = sparse_diag(d) * sparse_diff(sz,i,iscyclic[i])
  szt = collect(sz)
  
  if !iscyclic[i]
    # extx: extension operator
    szt[i] = szt[i]+1
    szt = (szt...)

    extx = sparse_trim(szt,i)'
    D = extx * D 
 else
    # shift back one grid cell along dimension i
    D = sparse_shift(sz,i,iscyclic[i])' * D
  end

  # add laplacian along dimension i  
  DD = DD + sparse_diff(szt,i,iscyclic[i]) * D  
end

ivol = .*(pmn...)

# Laplacian on regualar grid DD
DD = sparse_diag(ivol[:]) * DD


# Laplacian on grid with on sea points Lap
Lap = H * DD * H'


end

# Copyright (C) 2014,2016 Alexander Barth <a.barth@ulg.ac.be>
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
