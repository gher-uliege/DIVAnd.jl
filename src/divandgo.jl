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
  see divand_simple_example_go.jl
  
  Higher level with windowing etc
"""



function divandgo(mask,pmn,xi,x,f,Labs,epsilon2; otherargs...
                
                )

				
				
				
				
				
				
n=ndims(mask)
# Needed to make sure results are saved.
fi=0
s=0

# Need to check for cyclic boundaries

moddim=zeros(n);

kwargs_dict = Dict(otherargs)

if haskey(kwargs_dict,:moddim)
moddim=kwargs_dict[:moddim]
end

# DOES NOT YET WORK WITH PERIODIC DOMAINS OTHER THAN TO MAKE SURE THE DOMAIN IS NOT CUT 
# IN THIS DIRECTION. If adapation is done make sure the new moddim is passed to divandrun
# General approach in this future case prepare window indexes just taking any range including negative value and
# apply a mod(myindexes-1,size(mask)[i])+1 in direction i when extracting 
# for coordinates tuples of the grid (xin,yin, .. )  and data (x,y)
# in the direction, shift coordinates and apply modulo mod(x-x0+L/2,L)
#
#



# Analyse rations l/dx etc
Lpmnrange = divand_Lpmnrange(pmn,Labs)

# Create list of windows, steps for the coarsening during preconditioning and mask for lengthscales to decoupled directions during preconditioning
windowlist,csteps,lmask = divand_cutter(Lpmnrange,size(mask),moddim)


# For parallel version declare SharedArray(Float,size(mask)) instead of zeros() ? ? and add a @sync @parallel in front of the for loop ?
# Seems to work with an addprocs(2); @everywhere using divand to start the main program. To save space use Float32 ?
#fi=zeros(size(mask));
fi=SharedArray(Float64,size(mask));
@sync @parallel for iwin=1:size(windowlist)[1]
 
 iw1=windowlist[iwin][1]
 iw2=windowlist[iwin][2]
 isol1=windowlist[iwin][3]
 isol2=windowlist[iwin][4]
 istore1=windowlist[iwin][5]
 istore2=windowlist[iwin][6]

warn("Test window $iw1 $iw2 $isol1 $isol2 $istore1 $istore2 ")


windowpoints=([iw1[i]:iw2[i] for i in 1:n]...);



#################################################
# Need to check how to work with aditional constraints...
#################################################

#################################
# Search for velocity argument:
jfound=0
for j=1:size(otherargs)[1]
  if otherargs[j][1]==:velocity
  jfound=j
  break
  end
end


warn("There is an advection constraint; make sure the window sizes are large enough for the increased correlation length")



if jfound>0
# modify the parameter
   otherargsw=deepcopy(otherargs)
   otherargsw[jfound]=(:velocity,([ x[windowpoints...] for x in otherargs[jfound][2] ]...))
   else
   otherargsw=otherargs
end







# If C is square then maybe just take the sub-square corresponding to the part taken from x hoping the constraint is a local one ?
# 





# If C projects x on a low dimensional vector: maybe C'C x-C'd as a constraint, then pseudo inverse and woodbury to transform into a similar constraint but on each subdomain 
# Would for example replace a global average constraint to be replaced by the same constraint applied to each subdomain. Not exact but not too bad neither






fw=0
s=0
# Verify if a direct solver was requested from the demain decomposer
if sum(csteps)>0
fw,s=divandjog(mask[windowpoints...],([ x[windowpoints...] for x in pmn ]...),([ x[windowpoints...] for x in xi ]...),x,f,Labs,epsilon2,csteps,lmask; otherargsw...)
                else
fw,s=divandrun(mask[windowpoints...],([ x[windowpoints...] for x in pmn ]...),([ x[windowpoints...] for x in xi ]...),x,f,Labs,epsilon2; otherargsw...)				
end

# maskb=mask[windowpoints...]
# pmnb=([ x[windowpoints...] for x in pmn ]...)
# xib=([ x[windowpoints...] for x in xi ]...)

# fw,s=divandrun(maskb,pmnb,xib,x,f,Labs,epsilon2; otherargs...)
# Now error fields
# Cpme: just run and take out same window

# AEXERR: just run and take out same window


# EXERR: P only to points on the inner grid, not the overlapping one !


# End error fields



# Do similar things with other properties if asked for

# or better, replace the divandrun call by a divandrunsmart call with these parameters ? Would allow exploitation of optimised versions of the solver ?






# copy, deepcopy or just = ???



windowpointssol=([isol1[i]:isol2[i] for i in 1:n]...);
windowpointsstore=([istore1[i]:istore2[i] for i in 1:n]...);

#fi[istore1[1]:istore2[1],istore1[2]:istore2[2]]= fw[isol1[1]:isol2[1],isol1[2]:isol2[2]];

fi[windowpointsstore...]= fw[windowpointssol...];


end

# When finished apply an nd filtering to smooth possible edges, particularly in error fields.

# For the moment s is not defined ?
return fi,s 



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
