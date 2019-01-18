
function DIVAndrun(operatortype,mask::BitArray{N},pmnin,xiin,x,f::Vector{T},lin,epsilon2;
                   velocity = (),
                   primal::Bool = true,
                   factorize = true,
                   tol = 1e-6,
                   maxit = 100,
                   minit = 0,
                   constraints = (),
                   inversion = :chol,
                   moddim = [],
                   fracindex = Matrix{T}(undef,0,0),
                   alpha = [],
                   keepLanczosVectors = 0,
                   compPC = DIVAnd_pc_none,
                   progress = (iter,x,r,tol2,fun,b) -> nothing,
                   fi0 = zeros(size(mask)),
                   f0 = zeros(size(f)),
                   alphabc = 1.0,
                   scale_len = true,
                   btrunc=[],
				   MEMTOFIT=16.,
				   topographyforfluxes = (),
				   fluxes = (),
				   epsfluxes = 0,
				   RTIMESONESCALES=(),
				   QCMETHOD=()
                   ) where {N,T}

    # check pmn .* len > 4
    checkresolution(mask,pmnin,lin)

    pmn,xi,len = DIVAnd_bc_stretch(mask,pmnin,xiin,lin,moddim,alphabc)


    # observation error covariance (scaled)
    # Note: iB is scaled such that diag(inv(iB)) is 1 far from the
    # boundary
    # For testing this version of alphabc deactivate the other one
    s = DIVAnd_background(operatortype,mask,pmn,len,alpha,moddim,scale_len,[]; btrunc=btrunc);

    # check inputs
    if !any(mask[:])
        @warn "No sea points in mask, will return NaN";
        return fill(NaN,size(mask)),s
    end

    s.betap = 0;
    s.primal = primal;
    s.factorize = factorize;
    s.tol = tol;
    s.maxit = maxit;
    s.minit = minit;
    s.inversion = inversion;
    s.keepLanczosVectors = keepLanczosVectors;
    s.compPC = compPC;
    s.progress = progress

    #@info "Creating observation error covariance matrix"
    R = DIVAnd_obscovar(epsilon2,length(f));

    # add observation constrain to cost function
    #@info "Adding observation constraint to cost function"
    obscon = DIVAnd_obs(s,xi,x,f,R,fracindex)

    s = DIVAnd_addc(s,obscon);

    # add advection constraint to cost function
    if !isempty(velocity)
        #@info "Adding advection constraint to cost function"
        velcon = DIVAnd_constr_advec(s,velocity)
        s = DIVAnd_addc(s,velcon);
	end

	if !isempty(topographyforfluxes)
        #@info "Adding integral constraints"
        fluxcon = DIVAnd_constr_fluxes(s,topographyforfluxes,fluxes,epsfluxes,pmnin)
		s = DIVAnd_addc(s,fluxcon);
    end

    # add all additional constrains
    for i=1:length(constraints)
        #@info "Adding additional constrain - $(i)"
        s = DIVAnd_addc(s,constraints[i]);
    end

    # factorize a posteriori error covariance matrix
    # or compute preconditioner
    #@info "Factorizing a posteriori error covariance matrix"
    DIVAnd_factorize!(s);

    # @info "Solving..."
    fi0_pack = statevector_pack(s.sv,(fi0,))[:,1]

    #@code_warntype DIVAnd_solve!(s,fi0_pack,f0)
    fi = DIVAnd_solve!(s,fi0_pack,f0;btrunc=btrunc) :: Array{T,N}

    # @info "Done solving"
    return fi,s
end


function DIVAndrun(mask::Array{Bool,N},args...; kwargs...) where N
    return DIVAndrun(convert(BitArray{N},mask),args...; kwargs...)
end


# the same as DIVAndrun, but just return the field fi
DIVAndrunfi(args...) = DIVAndrun(args...)[1]



"""
    DIVAndrun(mask,pmn,xi,x,f,len,epsilon2; <keyword arguments>)

Perform an n-dimensional variational analysis of the observations `f` located at
the coordinates `x`. The array `fi` represent the interpolated field at the grid
defined by the coordinates `xi` and the scales factors `pmn`.

# Input:
* `mask`: binary mask delimiting the domain. true is inside and false outside.
For oceanographic application, this is the land-sea mask where sea is true and land is false.

* `pmn`: scale factor of the grid. pmn is a tuple with n elements. Every
   element represents the scale factor of the corresponding dimension. Its
   inverse is the local resolution of the grid in a particular dimension.
   For example, in two dimensions, `pmn` is a tuple `(pm,pn)` where `pm` is
   the inverse of the local resolution in first dimension and `pn` is the the inverse
   of the local resolution in second dimension.

*  `xi`: tuple with n elements. Every element represents a coordinate
  of the final grid on which the observations are interpolated.

* `x`: tuple with n elements. Every element represents a coordinate of
  the observations.

* `f`: value of the observations *minus* the background estimate (vector of
  `m` elements where `m` is the number of observations). See also note.

* `len`: tuple with n elements. Every element represents the correlation length for a given dimension.

* `epsilon2`: error variance of the observations (normalized by the error variance of the background field). `epsilon2` can be a scalar (all observations have the same error variance and their errors are decorrelated), a vector (all observations can have a different error variance and their errors are decorrelated) or a matrix (all observations can have a different error variance and their errors can be correlated). If `epsilon2` is a scalar, it is thus the *inverse of the signal-to-noise ratio*.

# Optional input arguments specified as keyword arguments

* `velocity`: velocity of advection constraint. The default is
       no-advection constraint

* `alpha`: alpha is vector of coefficients multiplying various terms in the
       cost function. The first element multiplies the norm.
       The other i-th element of alpha multiplies the (i+1)-th derivative.
       Per default, the highest derivative is m = ceil(1+neff/2) where neff is the
       effective dimension of the problem (the number of dimensions with a nonzero
       correlation length) and `ceil` is the ceiling function (rounding up).


       The values of alpha is the (m+1)th row of the Pascal triangle:
          m=0         1
          m=1       1   1
          m=1     1   2   1     (n=1,2)
          m=2   1   3   3   1   (n=3,4)
          ...

* `constraints`: a structure with user specified constraints (see `DIVAnd_addc`).

* `moddim`: modulo for cyclic dimension (vector with n elements).
     Zero is used for non-cyclic dimensions. One should not include a boundary
     zone (sometimes called a ghost zone or halo) for cyclic dimensions.
     For example if the first dimension
     is cyclic, then the grid point corresponding to `mask[1,j]` should be
     between `mask[end,1]` (left neighbor) and `mask[2,j]` (right neighbor).

* `fracindex`: fractional indices (n-by-m array). If this array is specified,
     then x and xi are not used.

* `inversion`: direct solver (:chol for Cholesky factorization) or an
     interative solver (:pcg for preconditioned conjugate gradient [1]) can be
     used.

* `compPC`: function that returns a preconditioner for the primal formulation
     if inversion is set to 'pcg'. The function has the following arguments:

           fun = compPC(iB,H,R)

    where iB is the inverse background error covariance, H the observation
    operator and R the error covariance of the observation. The function `compPC` returns the
    preconditioner `fun(x,fx)` computing fx = `M \\ x` (the inverse of M times x)
    where `M` is a positive defined symmetric matrix [1].
    Effectively, the system E⁻¹ A (E⁻¹)ᵀ (E x) = E⁻¹ b is solved for (E x) where E Eᵀ = M.
    Ideally, M should this be similar to A, so that E⁻¹ A (E⁻¹)ᵀ is close to the identity matrix.

* `fi0`: starting field for iterative primal algorithm (same size as `mask`).

* `f0`: starting field for iterative dual algorithm (same size as the observations `f`).

* `operatortype`: Val{:sparse} for using sparse matrices (default) or Val{:MatFun} or using functions
    to define the constrains.

* `scale_len`: true (default) if the correlation length-scale should be scaled
    such that the analytical
    kernel reaches 0.6019072301972346 (besselk(1.,1.)) at the same distance
    than in 2D. The kernel behaves thus similar to
    the default kernel in two dimensions (alpha = [1,2,1]).

* `alphabc` : numerical value defining how the last grid points are stretched outward.
   If `alphabc` is 1, the default value mimics an infinite domain.
   To have previous behaviour of finite domain use alphabc equal to `0`.

* `btrunc` : if provided defines where to truncate the calculation of the
    covariance matrix B. Only values up and including alpha[btrunc] will be
    calculated. If the iterative solution is calculated, the missing terms will
    be calculated on the fly during the conjugate gradient calulcations. Default value is none and full covariance calculation.

# Output:
*  `fi`: the analysed field
*  `s`: structure with an array `s.P` representing the analysed error covariance

# Note:
  If zero is not a valid first guess for your variable (as it is the case for
  e.g. ocean temperature), you have to subtract the first guess from the
  observations before calling DIVAnd and then add the first guess back in.

# Example:
  see DIVAnd_simple_example.jl

# References
[1]  https://en.wikipedia.org/w/index.php?title=Conjugate_gradient_method&oldid=761287292#The_preconditioned_conjugate_gradient_method
"""
function DIVAndrun(mask::BitArray,pmnin,xiin,x,f::Vector{T},lin,epsilon2;
                   operatortype = Val{:sparse}, kwargs...) where T

    return DIVAndrun(operatortype,mask,pmnin,xiin,x,f,lin,epsilon2; kwargs...)
end



# Copyright (C) 2008-2017 Alexander Barth <barth.alexander@gmail.com>
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

# LocalWords:  fi DIVAnd pmn len diag CovarParam vel ceil moddim fracdim
