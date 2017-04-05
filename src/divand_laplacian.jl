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

    DD = spzeros(prod(sz),prod(sz))

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

        D = oper_diag(operatortype,d) * oper_diff(operatortype,sz,i,iscyclic[i])

        if !iscyclic[i]
            # extx: extension operator
            szt = sz
            tmp_szt = collect(szt)
            tmp_szt[i] = tmp_szt[i]+1
            szt = (tmp_szt...)

            extx = oper_trim(operatortype,szt,i)'
            D = extx * D
            DD = DD + oper_diff(operatortype,szt,i,iscyclic[i]) * D
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
    Lap = H * DD * H'


end



#----------------------------------------------------------

for N = 1:6
@eval begin
function divand_laplacian_prepare{T}(mask::BitArray{$N},
                      pmn::NTuple{$N,Array{T,$N}},
                      nu::NTuple{$N,Array{T,$N}})
    const sz = size(mask)
    const ivol = .*(pmn...)

    const nus = ntuple(i -> zeros(sz),$N)::NTuple{$N,Array{T,$N}}

    # This heavily uses macros to generate fast code 
    # In e.g. 3 dimensions
    # (@nref $N tmp i) corresponds to tmp[i_1,i_2,i_3]
    # (@nref $N nu_i l->(l==j?i_l+1:i_l)  corresponds to nu_i[i_1+1,i_2,i_3] if j==1

    @nexprs $N j->begin
    tmp = nus[j]

    @nloops $N i k->(k == j ? (1:sz[k]-1) : (1:sz[k]))  begin
        nu_i = nu[j]
        # stagger nu
        # e.g. 0.5 * (nu[1][2:end,:] + nu[1][1:end-1,:])
        (@nref $N tmp i) = 0.5 * ((@nref $N nu_i l->(l==j?i_l+1:i_l)) + (@nref $N nu_i i) )
        # e.g. (pmn[2][2:end,:]+pmn[2][1:end-1,:]) ./ (pmn[1][2:end,:]+pmn[1][1:end-1,:])
        @nexprs $N m->begin
        pm_i = pmn[m]
        if (m .== j)
           (@nref $N tmp i) *= ((@nref $N pm_i l->(l==j?i_l+1:i_l)) + (@nref $N pm_i i) )
        else
           (@nref $N tmp i) /= ((@nref $N pm_i l->(l==j?i_l+1:i_l)) + (@nref $N pm_i i) )
        end

        if !(@nref $N mask i) || !(@nref $N mask l->(l==j?i_l+1:i_l))
            (@nref $N tmp i) = 0
        end
        end
    end
    end

    return ivol,nus
end

function divand_laplacian_apply{T}(ivol,nus,x::Array{T,$N})::Array{T,$N}
    sz = size(x)
    Lx = zeros(sz)

    @inbounds @nloops $N i d->1:sz[d]  begin
        (@nref $N Lx i) = 0
        @nexprs $N d1->begin      
            tmp2 = nus[d1]
            
            if i_d1 < sz[d1]
                (@nref $N Lx i) += (@nref $N tmp2 i) * ((@nref $N x d2->(d2==d1?i_d2+1:i_d2)) - (@nref $N x i))
            end
            
            if i_d1 > 1
                (@nref $N Lx i) -= (@nref $N tmp2 d2->(d2==d1?i_d2-1:i_d2)) * ((@nref $N x i) -  (@nref $N x d2->(d2==d1?i_d2-1:i_d2)))
            end
        end
        (@nref $N Lx i) *= (@nref $N ivol i)
    end
    
    return Lx
end

end # begin eval
end # for N = 1:6




# Copyright (C) 2014,2016,2017 Alexander Barth <a.barth@ulg.ac.be>
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
