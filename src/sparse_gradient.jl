for OT in [:sparse,:MatFun]
    ot = Val{Symbol(OT)}

    @eval begin
        @doc """
    Dx1,Dx2,...,Dxn = DIVAnd.DIVAnd_gradient(operatortype,mask,pmn,iscyclic)

Form the gradient using finite differences in all n-dimensions.
`mask` is a binary mask delimiting the domain. 1 is inside and 0 outside.
For oceanographic application, this is the land-sea mask.
`pmn` is a tuple of arrays with the scale factor of the grid.
The output `Dx1,Dx2,...,Dxn` are sparse matrices represeting a gradient along
different dimensions.
"""
        function DIVAnd_gradient(::Type{$ot},mask::AbstractArray{Bool,N},pmn::NTuple{N,AbstractArray{T,N}},iscyclic = falses(ndims(mask))) where {N,T}

            H = oper_pack($ot,mask)
            sz = size(mask)

            out = ntuple(i ->
                         begin
                         # staggering operator
                         S = oper_stagger($ot,sz,i,iscyclic[i])

                         # mask for staggered variable
                         m = (S * mask[:]) .== 1

                         d = m .* (S * pmn[i][:])

                         return oper_pack($ot,m) * oper_diag($ot,d) * oper_diff($ot,sz,i,iscyclic[i]) * H'
                         end,
                         Val(N))
            return out
        end
    end
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
