"""
Compute a variational analysis of arbitrarily located observations.

fi,s = divandgo(mask,pmn,xi,x,f,len,epsilon2,...);

Perform an n-dimensional variational analysis of the observations `f` located at
the coordinates `x`. The array `fi` represent the interpolated field at the grid
defined by the coordinates `xi` and the scales factors `pmn`.

# Input:
* `mask`: binary mask delimiting the domain. true is inside and false outside. For oceanographic application, this is the land-sea mask.

* `pmn`: scale factor of the grid. pmn is a tuple with n elements. Every
       element represents the scale factor of the corresponding dimension. Its
       inverse is the local resolution of the grid in a particular dimension.

*  `xi`: tuple with n elements. Every element represents a coordinate
  of the final grid on which the observations are interpolated

* `x`: tuple with n elements. Every element represents a coordinate of
  the observations

* `f`: value of the observations *minus* the background estimate (m-by-1 array).
    (see note)

* `len`: correlation length

* `epsilon2`: error variance of the observations (normalized by the error variance of the background field). `epsilon2` can be a scalar (all observations have the same error variance and their errors are decorrelated), a vector (all observations can have a difference error variance and their errors are decorrelated) or a matrix (all observations can have a difference error variance and their errors can be correlated). If `epsilon2` is a scalar, it is thus the *inverse of the signal-to-noise ratio*.

# Optional input arguments specified as keyword arguments

* `velocity`: velocity of advection constraint. The default is
       no-advection constraint

* `alpha`: alpha is vector of coefficients multiplying various terms in the
       cost function. The first element multiplies the norm.
       The other i-th element of alpha multiplies the (i+1)-th derivative.
       Per default, the highest derivative is m = ceil(1+n/2) where n is the
       dimension of the problem.

       The values of alpha is the (m+1)th row of the Pascal triangle:
          m=0         1
          m=1       1   1
          m=1     1   2   1     (n=1,2)
          m=2   1   3   3   1   (n=3,4)
          ...

* `diagnostics`: 0 or 1 turns diagnostic and debugging information on (1) or
       off (0, default). If on, they will be returned as the last output
       argument

* `EOF`, EOF: sub-space constraint. Orthogonal (EOF' WE^2 EOF = I) (units of
       EOF: m^(-n/2))

* `EOF_scaling`, EOF_scaling: (dimensional)

* `constraint`: a structure with user specified constrain

* `moddim`: modulo for cyclic dimension (vector with n elements).
     Zero is used for non-cyclic dimensions. Halo points should
     not be included for cyclic dimensions. For example if the first dimension
     is cyclic, then the grid point corresponding to mask(1,j) should be
     between mask(end,1) (left neighbor) and mask(2,j) (right neighbor)

* `fracdim`: fractional indices (n-by-m array). If this array is specified,
     then x and xi are not used.

* `inversion`: direct solver ('chol' for Cholesky factorization) or a
     interative solver ('pcg' for preconditioned conjugate gradient) can be
     used.

* `compPC`: function that returns a preconditioner for the primal formulation
     if inversion is set to 'pcg'. The function has the following arguments:

           [M1,M2] = compPC(iB,H,R)

    where iB is the inverse background error covariance, H the observation
    operator and R the error covariance of the observation. The used
    preconditioner M will be M = M1 * M2 (or just M = M1 if M2 is empty).
    Per default a modified incomplete Cholesky factorization will be used a
    preconditioner.

# Output:
*  `fi`: the analysed field
*  `s`: structure with an array `s.P` representing the analysed error covariance

# Note:
  If zero is not a valid first guess for your variable (as it is the case for
  e.g. ocean temperature), you have to subtract the first guess from the
  observations before calling divand and then add the first guess back in.

# Example:

  Higher level with recursive call with iterative method etc
"""



function divandjog(mask,pmn,xi,x,f,Labs,epsilon2; otherargs...
                
                )

				
				
				
				
				
				
n=ndims(mask)

# Need to check for cyclic boundaries

moddim=zeros(1,n);

kwargs_dict = Dict(otherargs)
@show itiscyclic=haskey(kwargs_dict, :moddim)

if itiscyclic
moddim=kwargs_dict[:moddim]
@show moddim  
end

@show moddim   

iscyclic = moddim .> 0

nsteps=divand_sampler(pmn,Labs);

@show nsteps
@show(sum(nsteps))

if sum(nsteps)>n 

####### use preconditionner methods




#######################################
# HOW TO TREAT COARSE GRID NOT COVERING FINE GRID ?
# Artificially add halo or also coordinates of last points in fine grid if the last step is beyond ?
#######################################

coarsegridpoints=([1:nsteps[i]:size(mask)[i] for i in 1:n]...);


xic=([ x[coarsegridpoints...] for x in xi ]...);

maskc=mask[coarsegridpoints...];

# Forget about land to start with
maskc=trues(size(maskc));

#Now scale pmn by the step factorsm if not possible during extraction this operation needs a copy
# Here hardcoded factor 3 
pmnc=([ (1.0/nsteps[i])*pmn[i][coarsegridpoints...] for i=1:length(pmn) ]...)

# Check if Labs is a tuple of tuple; in this case also subsample

	if isa(Labs,Tuple)
		if isa(Labs[1],Tuple)
		Labsc=([ x[coarsegridpoints...] for x in Labs ]...);   
           else
        Labsc=Labs;		   
		end
	else
	    Labsc=Labs;
	end



# Now prepare HI do go from the coarse grid to the fine grid. To do so
# interprete de fine grid coordinates as those of pseudo-obs and use divandtoos
#
# 
# Need the statevector strucure for the fine grid to go from grid to array
svf = statevector_init((mask,))

# For each coordinate in the tuplet xi, go from grid representation to tuplet
# to have the pseudo-data coordinates

#xfake=([statevector_pack(svf,(x,)) for x in xi]...)

#@show size(xfake[1])
# Create fractional indexes of these data points in the coarse grid
#Ic = localize_separable_grid(xfake,maskc,xic);



Ic = localize_separable_grid(([statevector_pack(svf,(x,)) for x in xi]...),maskc,xic);

#@show size(xfake[1])
# Create fractional indexes of these data points in the coarse grid

# Create HI
HI,outc,outbboxc = sparse_interp(maskc,Ic,iscyclic);

#@show maximum(xfake[1])

@show maximum(xic[1])


@time HI = HI * sparse_pack(maskc)';



####################################
# Need to look at constraints later
####################################






# Rune the coarse grid problem





#fw,s=divandrun(mask[windowpoints...],([ x[windowpoints...] for x in pmn ]...),([ x[windowpoints...] for x in xi ]...),x,f,Labs,epsilon2; otherargs...)





@time fc,sc=divandrun(maskc,pmnc,xic,x,f,Labsc,epsilon2; otherargs...)


@show Labsc

@show mean(fc)

# TEST OF Higher take the fc coarse solution, pack to to statevector form 
# using sc.sv
# Apply HI; this vector can also be used as a first guess for the PC

   xguess=HI*statevector_pack(sc.sv,(fc,));
   
@show mean(xguess)

# which you unpack using the statevector form of the fine grid for the TEST only

#   figuess,=statevector_unpack(svf,xguess)
# Do not know why I need to squeeze here if I want to return a gridded approximation
#   figuess=squeeze(figuess,ndims(figuess))
# Recover sc.P and define the conditionner

# tolerance on the gradient A x - b
tol = 5e-2


maxiter=10*Int(ceil(sqrt(size(HI)[1])))
@show maxiter

pcargs = [(:tol, tol),(:maxit,maxiter)]

@show size(HI)
#@show sc.P.factors[:PtL]*ones(size(HI'))

#jmPHI=sc.P.factors[:PtL]\copy(HI');

diagshift=0.01*(sqrt(size(HI)[1]/size(HI)[2])-1);
@show diagshift

function compPC(iB,H,R)
        return x -> diagshift*x+HI*(sc.P*(HI'*x));
	#     return jmPHI'*(jmPHI*x);
	#   return x->x;
end
# First guess is the HI* coarse solution

# HI*(sc.P*(HI'  *x ))  should be a good operator for the M-1 x operation in preconditionner ?
# Why do I need to take sc.P\ ??? So better use components of P*HI' ?


# Then run with normal resolution and preconditionner
fi,si=divandrun(mask,pmn,xi,x,f,Labs,epsilon2; otherargs...,pcargs...,inversion=:pcg,compPC = compPC, fi0 =xguess)

errfield=diagMtCM(sc.P,HI')

# 
#erri,=statevector_unpack(si,errfield)
# For error field based on coarse one, use divand_filter3 with ntimes=Int(ceil(mean(nsteps)))


# First test
return fi,si

#fc,sc,figuess,fi,si,errfield
#####################################################
# end of iterative version
#####################################################
else
#####################################################
# Run normal version
#####################################################

finter,sinter=divandrun(mask,pmn,xi,x,f,Labs,epsilon2; otherargs...)

return finter,sinter
end
##########################

end

# Copyright (C) 2008-2017 Alexander Barth <barth.alexander@gmail.com>
#                         Jean-Marie Beckers   <JM.Beckers@ulg.ac.be>
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

# LocalWords:  fi divand pmn len diag CovarParam vel ceil moddim fracdim
