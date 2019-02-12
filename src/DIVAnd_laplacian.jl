"""
 Create the laplacian operator.

 Lap = DIVAnd_laplacian(mask,pmn,nu,iscyclic)

 Form a Laplacian using finite differences
 assumes that gradient is zero at "coastline"

 Input:
   mask: binary mask delimiting the domain. 1 is inside and 0 outside.
         For oceanographic application, this is the land-sea mask.
   pmn: scale factor of the grid.
   nu: diffusion coefficient of the Laplacian
      field of the size mask or cell arrays of fields

 Output:
   Lap: sparce matrix represeting a Laplacian

"""
function DIVAnd_laplacian(operatortype,mask::AbstractArray{Bool,N},pmn,nu::Float64,iscyclic) where N
    nu_ = ntuple(fill(nu,size(mask)),Val(N))

    return DIVAnd_laplacian(operatortype,mask,pmn,nu_,iscyclic)

end

function DIVAnd_laplacian(operatortype,mask,pmn,nu::Array{Float64,N},iscyclic) where N
    nu_ = ntuple(i -> nu[i],Val(N))
    return DIVAnd_laplacian(operatortype,mask,pmn,nu_,iscyclic)
end

function DIVAnd_laplacian(operatortype,mask,pmn,nu::Tuple{Vararg{Number,N}},iscyclic) where N
    nu_ = ntuple(i -> fill(nu[i],size(mask)),Val(N))
    return DIVAnd_laplacian(operatortype,mask,pmn,nu_,iscyclic)
end

function DIVAnd_laplacian(operatortype,mask,pmn,nu::Tuple{Vararg{AbstractArray{T,n},n}},iscyclic;
                          coeff_laplacian::Vector{Float64} = ones(ndims(mask))
                          ) where {T,n}
    sz = size(mask)

    # extraction operator of sea points
    H = oper_pack(operatortype,mask)

	# Already include final operation *H' on DD to reduce problem size so resized the matrix
    # DD = spzeros(prod(sz),prod(sz))
	DD = spzeros(size(H,2),size(H,1))

    for i=1:n
        if coeff_laplacian[i] != 0
            # operator for staggering in dimension i
            S = oper_stagger(operatortype,sz,i,iscyclic[i])

            # d = 1 for interfaces surounded by water, zero otherwise
            d = zeros(Float64,size(S,1))
            d[(S * mask[:]) .== 1] .= 1.

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

            d = coeff_laplacian[i] * (d .* (S * nu[i][:]))
            #d = d .* (S * nu[i][:])

            # Flux operators D
            # zero at "coastline"

            #D = oper_diag(operatortype,d) * oper_diff(operatortype,sz,i,iscyclic[i])
            # Already include final operation to reduce problem size in real problems with land mask

            D = oper_diag(operatortype,d) * (oper_diff(operatortype,sz,i,iscyclic[i])*H')

            # add laplacian along dimension i
            if !iscyclic[i]
                # extx: extension operator
                tmp_szt = collect(sz)
                tmp_szt[i] = tmp_szt[i]+1
                szt = NTuple{n,Int}(tmp_szt)
                #szt = ntuple(j -> (j == i ? size(mask,j)+1 : size(mask,j) ),Val(n))

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
        end
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

    return Lap
end


#----------------------------------------------------------

for N = 1:6
@eval begin
function DIVAnd_laplacian_prepare(mask::BitArray{$N},
                      pmn::NTuple{$N,Array{T,$N}},
                      nu::NTuple{$N,Array{T,$N}}) where T
    sz = size(mask)
    ivol = .*(pmn...)

    nus = ntuple(i -> zeros(T,sz),$N)::NTuple{$N,Array{T,$N}}

    # This heavily uses macros to generate fast code
    # In e.g. 3 dimensions
    # (@nref $N tmp i) corresponds to tmp[i_1,i_2,i_3]
    # (@nref $N nu_i l->(l==j?i_l+1:i_l)  corresponds to nu_i[i_1+1,i_2,i_3] if j==1

    # loop over all dimensions to create
    # nus[1] (nu stagger in the 1st dimension)
    # nus[2] (nu stagger in the 2nd dimension)
    # ...

    @nexprs $N j->begin
    tmp = nus[j]

    # loop over all spatio-temporal dimensions
    @nloops $N i k->(k == j ? (1:sz[k]-1) : (1:sz[k]))  begin
        nu_i = nu[j]
        # stagger nu
        # e.g. 0.5 * (nu[1][2:end,:] + nu[1][1:end-1,:])
        (@nref $N tmp i) = 0.5 * ((@nref $N nu_i l->(l==j ? i_l+1 : i_l)) + (@nref $N nu_i i) )
        # e.g. (pmn[2][2:end,:]+pmn[2][1:end-1,:]) ./ (pmn[1][2:end,:]+pmn[1][1:end-1,:])
        @nexprs $N m->begin
        pm_i = pmn[m]
        if (m .== j)
            (@nref $N tmp i) *= 0.5 * ((@nref $N pm_i l->(l==j ? i_l+1 : i_l)) + (@nref $N pm_i i) )
        else
            (@nref $N tmp i) /= 0.5 * ((@nref $N pm_i l->(l==j ? i_l+1 : i_l)) + (@nref $N pm_i i) )
        end

        if !(@nref $N mask i) || !(@nref $N mask l->(l==j ? i_l+1 : i_l))
            (@nref $N tmp i) = 0
        end
        end
    end
    end

    return ivol,nus
end


function DIVAnd_laplacian_apply!(ivol,nus,x::AbstractArray{T,$N},Lx::AbstractArray{T,$N}) where T
    sz = size(x)
    Lx[:] .= 0

    @inbounds @nloops $N i d->1:sz[d]  begin
        (@nref $N Lx i) = 0
        @nexprs $N d1->begin
            tmp2 = nus[d1]

            if i_d1 < sz[d1]
                @inbounds (@nref $N Lx i) += (@nref $N tmp2 i) * (
                    (@nref $N x d2->(d2==d1 ? i_d2+1 : i_d2)) - (@nref $N x i))
            end

            if i_d1 > 1
                @inbounds (@nref $N Lx i) -= (@nref $N tmp2 d2->(d2==d1 ? i_d2-1 : i_d2)) * (
                    (@nref $N x i) -  (@nref $N x d2->(d2==d1 ? i_d2-1 : i_d2)))
            end
        end
        (@nref $N Lx i) *= (@nref $N ivol i)
    end
end


function DIVAnd_laplacian_apply(ivol,nus,x::AbstractArray{T,$N})::AbstractArray{T,$N} where T
    Lx = similar(x)
    DIVAnd_laplacian_apply!(ivol,nus,x,Lx)
    return Lx
end



end # begin eval
end # for N = 1:6




function _derivative2n!(dim,mask,pm,len,va,D,Rpre,n,Rpost)
    for Ipost in Rpost
        for i = 2:n-1
            for Ipre in Rpre
                if mask[Ipre,i-1,Ipost] && mask[Ipre,i,Ipost] && mask[Ipre,i+1,Ipost]
                    D[Ipre,i,Ipost] += len[Ipre,i,Ipost]^2 * pm[Ipre,i,Ipost]^2 * (va[Ipre,i-1,Ipost] - 2*va[Ipre,i,Ipost] + va[Ipre,i+1,Ipost])
                end
            end
        end
    end
end

function derivative2n!(dim::Integer,mask,pmn,len,va,D)
    Rpre = CartesianIndices(size(va)[1:dim-1])
    Rpost = CartesianIndices(size(va)[dim+1:end])

    _derivative2n!(dim,mask,pmn[dim],len[dim],va,D,Rpre,size(va,dim),Rpost)
    return D
end

derivative2n(dim,mask,pmn,len,va) = derivative2n!(dim,mask,pmn,len,va,zeros(size(mask)))



function _sparse_derivative2n!(dim,mask,pm,len,Rpre,n,Rpost)
    S = sparse([], [], Float64[],length(mask),length(mask))
    linindex = LinearIndices(mask)

    for Ipost in Rpost
        for i = 2:n-1
            for Ipre in Rpre
                if mask[Ipre,i-1,Ipost] && mask[Ipre,i,Ipost] && mask[Ipre,i+1,Ipost]
                    coeff = len[Ipre,i,Ipost]^2 * pm[Ipre,i,Ipost]^2

                    S[linindex[Ipre,i,Ipost],linindex[Ipre,i-1,Ipost]] += coeff
                    S[linindex[Ipre,i,Ipost],linindex[Ipre,i  ,Ipost]] += -2*coeff
                    S[linindex[Ipre,i,Ipost],linindex[Ipre,i+1,Ipost]] += coeff
                end
            end
        end
    end
    return S
end

function sparse_derivative2n(dim::Integer,mask,pmn,len)
    Rpre = CartesianIndices(size(mask)[1:dim-1])
    Rpost = CartesianIndices(size(mask)[dim+1:end])

    S = _sparse_derivative2n!(dim,mask,pmn[dim],len[dim],Rpre,size(mask,dim),Rpost)
    return S
end



# Copyright (C) 2014,2016-2019 Alexander Barth 		<a.barth@ulg.ac.be>
#                              Jean-Marie Beckers 	<jm.beckers@ulg.ac.be>
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
