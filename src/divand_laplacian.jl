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
#   Lap: sparce matrix representing a Laplacian
#
#

function divand_laplacian(operatortype,mask,pmn,nu::Float64,iscyclic)
    n = ndims(mask)
    nu_ = ((fill(nu,size(mask)) for i = 1:n)...)

    return divand_laplacian(operatortype,mask,pmn,nu_,iscyclic)

end

function divand_laplacian{n}(operatortype,mask,pmn,nu::Array{Float64,n},iscyclic)
    nu_ = ((nu for i = 1:n)...)
    return divand_laplacian(operatortype,mask,pmn,nu_,iscyclic)
end

function divand_laplacian{n}(operatortype,mask,pmn,nu::Tuple{Vararg{Number,n}},iscyclic)
    nu_ = ((fill(nui,size(mask)) for nui = nu)...)
    return divand_laplacian(operatortype,mask,pmn,nu_,iscyclic)
end

function divand_laplacian{n}(operatortype,mask,pmn,nu::Tuple{Vararg{Any,n}},iscyclic)

    sz = size(mask)

    # default float type

    # extraction operator of sea points
    H = oper_pack(operatortype,mask)
	
    sz = size(mask)
	# Already include final operation *H' on DD to reduce problem size so resized the matrix
    # DD = spzeros(prod(sz),prod(sz))
	DD = spzeros(size(H)[2],size(H)[1])

    for i=1:n
        # operator for staggering in dimension i
        S = oper_stagger(operatortype,sz,i,iscyclic[i])

        # d = 1 for interfaces surounded by water, zero otherwise
        d = zeros(Float64,size(S,1))
        d[(S * mask[:]) .== 1] = 1.

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

        #D = oper_diag(operatortype,d) * oper_diff(operatortype,sz,i,iscyclic[i])
		# Already include final operation to reduce problem size in real problems with land mask
        D = oper_diag(operatortype,d) * (oper_diff(operatortype,sz,i,iscyclic[i])*H')
		
        if !iscyclic[i]
            # extx: extension operator
            szt = sz
            tmp_szt = collect(szt)
            tmp_szt[i] = tmp_szt[i]+1
            szt = (tmp_szt...)

            extx = oper_trim(operatortype,szt,i)'
			# Tried to save an explicitely created intermediate matrix D
            #D = extx * D
            #DD = DD + oper_diff(operatortype,szt,i,iscyclic[i]) * D
			#@show size(extx),size(D),typeof(extx),typeof(D)
			#@show extx
			DD = DD + oper_diff(operatortype,szt,i,iscyclic[i]) * (extx * D)
        else
            # shift back one grid cell along dimension i
            D = oper_shift(operatortype,sz,i,iscyclic[i])' * D
            DD = DD + oper_diff(operatortype,sz,i,iscyclic[i]) * D
        end

        # add laplacian along dimension i
    end

    ivol = .*(pmn...)

    # Laplacian on regualar grid DD
    DD = oper_diag(operatortype,ivol[:]) * DD


    # Laplacian on grid with on sea points Lap
    #Lap = H * DD * H'
	# *H' already done before
	Lap = H * DD
#@show Lap
# Possible test to force symmetric Laplacian
#    Lap=0.5*(Lap+Lap')

end

# Copyright (C) 2014,2017 Alexander Barth 		<a.barth@ulg.ac.be>
#                         Jean-Marie Beckers 	<jm.beckers@ulg.ac.be>
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
