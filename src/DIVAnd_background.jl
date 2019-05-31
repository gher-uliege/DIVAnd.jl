"""
Form the inverse of the background error covariance matrix.
s = DIVAnd_background(mask,pmn,Labs,alpha,moddim)
Form the inverse of the background error covariance matrix with
finite-difference operators on a curvilinear grid
# Input:
* mask: binary mask delimiting the domain. 1 is inside and 0 outside.
        For oceanographic application, this is the land-sea mask.
* pmn: scale factor of the grid.
* Labs: correlation length
* alpha: a dimensional coefficients for norm, gradient, laplacian,...
     alpha is usually [1,2,1] in 2 dimensions.
# Output:
*  s: stucture containing
    * s.iB: inverse of the background error covariance
    * s.L: spatial average correlation length
    * s.n: number of dimenions
    * s.coeff: scaling coefficient such that the background variance diag(inv(iB)) is one far away from the boundary.
"""
function DIVAnd_background(operatortype,mask,pmn,Labs,alpha,moddim,scale_len = true,mapindex = [];
                           btrunc = [],
                           coeff_laplacian::Vector{Float64} = ones(ndims(mask)),
                           coeff_derivative2::Vector{Float64} = zeros(ndims(mask))
                           )

    # number of dimensions
    n = ndims(mask)

    Labs = len_harmonize(Labs,mask)

    neff, alpha = alpha_default(Labs,alpha)

    sz = size(mask)

    if isempty(moddim)
        moddim = zeros(n)
    end

    iscyclic = moddim .> 0

    # scale iB such that the diagonal of inv(iB) is 1 far from
    # the boundary
    # we use the effective dimension neff to take into account that the
    # correlation length-scale might be zero in some directions


    coeff = 1
    len_scale = 1
    try
        coeff,K,len_scale = DIVAnd_kernel(neff,alpha)
    catch err
        if isa(err, DomainError)
            @warn "no scaling for alpha=$(alpha)"
        else
            rethrow(err)
        end
    end

    if scale_len
        # scale Labs by len_scale so that all kernels are similar
        Labs = ntuple(i -> Labs[i]/len_scale,n)
    end


    # mean correlation length in every dimension
    Ld = [mean(L) for L in Labs]
    neff = sum(Ld .> 0)

    @debug "effective number of dimensions (neff): $neff"

    # geometric mean
    geomean(v) = prod(v)^(1/length(v))
    L = geomean(Ld[Ld .> 0])

    alphabc = 0

    s,D = DIVAnd_operators(operatortype,mask,pmn,([L.^2 for L in Labs]...,),
                           iscyclic,mapindex,Labs;
                           coeff_laplacian = coeff_laplacian,
                           )

    # D is laplacian (a dimensional, since nu = Labs.^2)
    sv = s.sv
    n = s.n


    # Labsp: 1st index represents the dimensions
    #Labsp = permute(Labs,[n+1 1:n])
    #pmnp = permute(pmn,[n+1 1:n])

    # mean correlation length in every dimension

    # # geometric mean
    # geomean(v) = prod(v)^(1/length(v))
    # L = geomean(Ld[Ld .> 0])


    # norm taking only dimension into account with non-zero correlation
    # WE: units length^(neff/2)

    d = .*(pmn[Ld .> 0]...)

	WE = oper_diag(operatortype,statevector_pack(sv,(1 ./ sqrt.(d),))[:,1])



	Ln = prod(Ld[Ld .> 0])



#if any(Ld <= 0)
#   pmnd = mean(reshape(pmnp,[n sv.numels_all]),2)
#   #Ln = Ln * prod(pmnd(Ld <= 0))
#end

	coeff = coeff * Ln # units length^n

	pmnv = hcat([pm[:] for pm in pmn]...)

	pmnv[:,findall(Ld .== 0)] .= 1

# staggered version of norm

	for i=1:n
		S = sparse_stagger(sz,i,iscyclic[i])
		ma = (S * mask[:]) .== 1
        d = @static if VERSION >= v"0.7.0-beta.0"
		    sparse_pack(ma) * (prod(S * pmnv,dims=2)[:,1])
        else
		    sparse_pack(ma) * (prod(S * pmnv,2)[:,1])
        end
		d = 1 ./ d
		s.WEs[i] = oper_diag(operatortype,sqrt.(d))
	end

# staggered version of norm scaled by length-scale

#s.Dxs = []
	for i=1:n
		Li2 = Labs[i][:].^2

		S = sparse_stagger(sz,i,iscyclic[i])

		# mask for staggered variable
		m = (S * mask[:]) .== 1

		tmp = sparse_pack(m) * sqrt.(S*Li2[:])
		s.WEss[i] = oper_diag(operatortype,tmp) * s.WEs[i]
		#  s.Dxs[i] = sparse_diag(sqrt(tmp)) * s.Dx[i]
	end

# adjust weight of halo points
	if !isempty(mapindex)
		# ignore halo points at the center of the cell

		WE = oper_diag(operatortype,s.isinterior) * WE

    # divide weight be two at the edged of halo-interior cell
    # weight of the grid points between halo and interior points
    # are 1/2 (as there are two) and interior points are 1
		for i=1:n
			s.WEs[i] = oper_diag(operatortype,sqrt.(s.isinterior_stag[i])) * s.WEs[i]
		end
	end

	s.WE = WE
	s.coeff = coeff
	# number of dimensions
	s.n = n

	# mean correlation legth
	s.Ld = Ld
    s.coeff_derivative2 = coeff_derivative2
    s.coeff_laplacian = coeff_laplacian

	iB = DIVAnd_background_components(
        s,D,alpha,btrunc=btrunc,coeff_derivative2=coeff_derivative2)

    # second order derivative background constraint without cross-terms
    pack = sparse_pack(mask)
    for i = 1:n
        if coeff_derivative2[i] != 0.
            S = sqrt(coeff_derivative2[i]) * s.WE * pack * DIVAnd.sparse_derivative2n(i,mask,pmn,Labs) * pack'
            iB += S' * S
        end
    end

	# inverse of background covariance matrix
	s.iB = iB


	#s.Ln = Ln

	s.moddim = moddim
	s.iscyclic = iscyclic

	s.alpha = alpha
	s.neff = neff
	s.WE = WE # units length^(n/2)

	return s
end

# Copyright (C) 2014, 2017 Alexander Barth 		<a.barth@ulg.ac.be>
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
