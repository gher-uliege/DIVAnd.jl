"""
Form the inverse of the background error covariance matrix.

s = divand_background(mask,pmn,Labs,alpha,moddim)

Form the inverse of the background error covariance matrix with
finite-difference operators on a curvilinear grid

# Input:
* mask: binary mask delimiting the domain. 1 is inside and 0 outside.
        For oceanographic application, this is the land-sea mask.

* pmn: scale factor of the grid.

* Labs: correlation length

* alpha: a dimensional coefficients for norm, gradient, laplacian,...
     alpha is usually [1 2 1] in 2 dimensions.

# Output:
*  s: stucture containing
    * s.iB: inverse of the background error covariance
    * s.L: spatial average correlation length
    * s.n: number of dimenions
    * s.coeff: scaling coefficient such that the background variance diag(inv(iB)) is one far away from the boundary.
"""
function divand_background(mask,pmn,Labs,alpha,moddim,mapindex = [];alphabc=2)

# number of dimensions
n = ndims(mask)
sz = size(mask)

if isempty(moddim)
  moddim = zeros(1,n)
end

iscyclic = moddim .> 0


if isa(Labs,Number)
    Labs = ((Labs * ones(size(mask)) for i=1:n)...)
elseif isa(Labs,Tuple)

    if isa(Labs[1],Number)
        Labs = ([Labs[i] * ones(size(mask)) for i = 1:n]...)
    end

    for i=1:n
        if !isequal(size(mask),size(Labs[i]))
            error("mask (%s) and correlation length (%s) have incompatible size",
                  formatsize(size(mask)),formatsize(size(Labs[i])))
        end
    end
end

if isempty(alpha)
  # kernel should has be continuous derivative

  # highest derivative in cost function
  m = Int(ceil(1+n/2))

  # alpha is the (m+1)th row of the Pascal triangle:
  # m=0         1
  # m=1       1   1
  # m=1     1   2   1
  # m=2   1   3   3   1
  # ...

  alpha = [binomial(m,k) for k = 0:m]
end

#if ~isequal([size(mask) n],size(pmn))
#  error('mask (#s) and metric (#s) have incompatible size',formatsize(size(mask)),formatsize(size(pmn)))
#end

s = divand_operators(mask,pmn,([_.^2 for _ in Labs]...),iscyclic,mapindex)
D = s.D # laplacian (a dimensional, since nu = Labs.^2)
sv = s.sv
n = s.n


# Labsp: 1st index represents the dimensions
#Labsp = permute(Labs,[n+1 1:n])
#pmnp = permute(pmn,[n+1 1:n])

# must handle the case when Labs is zero in some dimension
# thus reducing the effective dimension

# mean correlation length in every dimension
Ld = [mean(_) for _ in Labs]
neff = sum(Ld .> 0)
#L = mean(Ld(Ld > 0))

# geometric mean
geomean(v) = prod(v)^(1/length(v))
L = geomean(Ld[Ld .> 0])


# norm taking only dimension into account with non-zero correlation
# WE: units length^(neff/2)
#d = prod(pmnp(find(Ld > 0),:),1)'

d = .*(pmn[Ld .> 0]...)

#JMB This is probably the place to increase weight on anomaly constrained near borders.
# d is 1/volume of each grid box and an n dimensional array ?
# Use of CartesianRange to get boundaries of the n-dimensional box ?
# Or slicedime exploiting use by reference ??? in which the matrix from where the slice is "taken"
# is also updated ? Check in simplecase Nope does not work as reshape ... 
# view works ! view(A,2,1)[1,1]=1000

# For the moment hardcoded in 2D AND with average length scale instead of local one

# d=d./(alphabc.*pmn.*l)
# Dimension 1
#alphabc=0
@show alphabc
if alphabc>0

if n==1

d[1]=d[1]./(alphabc.*pmn[1][1].*Ld[1])
d[end]=d[end]./(alphabc.*pmn[1][end].*Ld[1])

end


if n==2

d[1,:]=d[1,:]./(alphabc.*pmn[1][1,:].*Ld[1])
d[end,:]=d[end,:]./(alphabc.*pmn[1][end,:].*Ld[1])

d[:,1]=d[:,1]./(alphabc.*pmn[2][:,1].*Ld[2])
d[:,end]=d[:,end]./(alphabc.*pmn[2][:,end].*Ld[2])

end



end

#/JMB


WE = sparse_diag(statevector_pack(sv,(1./sqrt.(d),))[:,1])


# scale iB such that the diagonal of inv(iB) is 1 far from
# the boundary
# we use the effective dimension neff to take into account that the
# correlation length-scale might be zero in some directions

Ln = prod(Ld[Ld .> 0])

#if any(Ld <= 0)
#   pmnd = mean(reshape(pmnp,[n sv.numels_all]),2)
#   #Ln = Ln * prod(pmnd(Ld <= 0))
#end


coeff = 1
try
    coeff,K = divand_kernel(neff,alpha)
    coeff = coeff * Ln # units length^n
catch err
    if isa(err, DomainError)
        warn("no scaling for alpha=$(alpha)")
    else
       rethrow(err)
    end
end

pmnv = cat(2,[_[:] for _ in pmn]...)

#JMB
# change here for the metric? check what cat does exactly then just modify pmnv !


pmnv[:,find(Ld == 0)] = 1

# staggered version of norm

for i=1:n
    S = sparse_stagger(sz,i,iscyclic[i])

    ma = (S * mask[:]) .== 1
    d = sparse_pack(ma) * (prod(S * pmnv,2)[:,1])
    d = 1./d
  s.WEs[i] = sparse_diag(sqrt.(d))
end

# staggered version of norm scaled by length-scale

#s.Dxs = []
for i=1:n
    Li2 = Labs[i][:].^2

    S = sparse_stagger(sz,i,iscyclic[i])

    # mask for staggered variable
    m = (S * mask[:]) .== 1

    tmp = sparse_pack(m) * sqrt.(S*Li2[:])
    s.WEss[i] = sparse_diag(tmp) * s.WEs[i]
    #  s.Dxs[i] = sparse_diag(sqrt(tmp)) * s.Dx[i]
end

# adjust weight of halo points
if !isempty(mapindex)
  # ignore halo points at the center of the cell

  WE = sparse_diag(s.isinterior) * WE

  # divide weight be two at the edged of halo-interior cell
    # weight of the grid points between halo and interior points
    # are 1/2 (as there are two) and interior points are 1
  for i=1:n
    s.WEs[i] = sparse_diag(sqrt.(s.isinterior_stag[i])) * s.WEs[i]
  end
end

s.WE = WE
s.coeff = coeff
# number of dimensions
s.n = n

iB_,iB = divand_background_components(s,alpha)

# inverse of background covariance matrix
s.iB = iB

# individual terms of background covariance matrix (corresponding alpha)
s.iB_ = iB_

#s.L = L
#s.Ln = Ln

s.moddim = moddim
s.iscyclic = iscyclic

s.alpha = alpha
s.neff = neff
s.WE = WE # units length^(n/2)

return s
end

# Copyright (C) 2014, 2016 Alexander Barth <a.barth@ulg.ac.be>
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
