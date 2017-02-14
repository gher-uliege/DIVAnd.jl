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

# Need to check for cyclic boundaries

moddim=zeros(n);

kwargs_dict = Dict(otherargs)
@show itiscyclic=haskey(kwargs_dict, :moddim)

if itiscyclic
moddim=kwargs_dict[:moddim]
@show moddim  
end

@show moddim   
# DOES NOT YET WORK WITH PERIODIC DOMAINS OTHER THAN TO MAKE SURE THE DOMAIN IS NOT CUT 
# IN THIS DIRECTION. If adapation is done make sure the new moddim is passed to divandrun
# General approach in this future case prepare window indexes just taking any range including negative value and
# apply a mod(myindexes-1,size(mask)[i])+1 in direction i when extracting 
# for coordinates tuples of the grid (xin,yin, .. )  and data (x,y)
# in the direction, shift coordinates and apply modulo mod(x-x0+L/2,L)
#
#

# Also there is a huge overhead in the test_divandgo case. Need to analyze


# Hardwired

# How wide is the overlap in terms of number of length scales
factoroverlap=2

# If possible use windows as large as possible but small enought to not reach the following problem size
# This needs adaptions to the dimensionality of the problem; the more dimensions the more expensive
# the problem is because of the bandwidth of the matrix. So typically

if n<3
biggestproblem=300*300 
end
if n==3
biggestproblem=50*50*50
end
if n==4
biggestproblem=100*100*3*12
end

# But increase it again if divandjog is called and can use subsampling...
laterscales=divand_sampler(pmn,Labs)
@show biggestproblem
biggestproblem=biggestproblem*prod(laterscales)
@show biggestproblem

# Fraction of domain size at which windowing is attempted if len is smaller
lfactor=0.2

pmnw=[];
maskw=[];
xw=[];
fw=[];
ndlast=0;

fi=zeros(size(mask));
# check inputs



Lscales=ones(n);
stepsize=ones(Int,n);
overlapping=ones(Int,n);
Lscalespmnmax=ones(n);

if isa(Labs,Number)
    Lscales=Labs*Lscales;
	for i=1:n
	 Lscalespmnmax[i]=Lscales[i]*maximum(pmn[i]);
	end
	
		
elseif isa(Labs,Tuple)

    if isa(Labs[1],Number)
        Lscales = [[Labs[i]  for i = 1:n]...]
		for i=1:n
	    Lscalespmnmax[i]=Lscales[i]*maximum(pmn[i]);
	    end
	else
     	for i=1:n
	    Lscalespmnmax[i]=maximum(Labs[i].*pmn[i]);
	    end
    end

end

@show Lscales
@show Lscalespmnmax

problemsize=1;
nwd=0
for i=1:n
# For each dimension try to define windows
# Assuming uniform metrics
# Default is no window:
stepsize[i]=size(mask)[i];
overlapping[i]=0;

# if length scale is small compared to domain size
if Lscalespmnmax[i]<   lfactor*size(mask)[i]
 if moddim[i]==0
 overlapping[i]=Int(ceil( factoroverlap*Lscalespmnmax[i]   ))
 problemsize=problemsize*overlapping[i]
 # Forced small window on zero
  if i==3
  overlapping[i]=1
  end
 
  nwd=nwd+1 
                              else
 problemsize=problemsize*size(mask)[i]		
 end 
else
problemsize=problemsize*size(mask)[i]		
end

# 
end

@show problemsize
@ show nwd
@ show overlapping

if nwd>0
  epsilon=(float(biggestproblem)/float(problemsize))^(1.0/nwd)-2.0
end
if epsilon<0
warn("SO what $epsilon $problemsize $nwd $overlapping")
epsilon=1E-6
end

for i=1:n
# if length scale is small compared to domain size
if Lscalespmnmax[i]<   lfactor*size(mask)[i]
if moddim[i]==0
stepsize[i]=Int(ceil( epsilon*factoroverlap*Lscalespmnmax[i]   ))
end
end
  if i==3
  stepsize[i]=1
  end
# 
end



# Depending on problem size increase or decrease stepsize for the windowed part






# Keep pointers to place where solution is found in output grid
# window corners in main problem
iw1=zeros(Int,n)
iw2=zeros(Int,n)

# Solution indexes in submodel
isol1=zeros(Int,n)
isol2=zeros(Int,n)
ij=ones(Int,n)

# Range for final storage
istore1=zeros(Int,n)
istore2=zeros(Int,n)


# Now test hardcoded for maximum 4 dimensions
ij=zeros(Int,4);
ssize=ones(Int,4);
lastp=ones(Int,4);

ssize[1:n]=stepsize[1:n];
lastp[1:n]=collect(size(mask));



sz = size(mask)
subsz = ([ceil(Int,sz[i] / stepsize[i]) for i = 1:ndims(mask)]...)
for cr in CartesianRange(subsz)
    ij = [[(cr[i]-1)*stepsize[i]+1  for i = 1:ndims(mask)]...]
        @show ij
# end 

# for il4=1:ssize[4]:lastp[4]
# for il3=1:ssize[3]:lastp[3]
# for il2=1:ssize[2]:lastp[2]
# for il1=1:ssize[1]:lastp[1]

# ij[1]=il1;
# ij[2]=il2;
# ij[3]=il3;
# ij[4]=il4;




# Generic code again
for nd=1:n
iw1[nd]=ij[nd]-overlapping[nd]
isol1[nd]=overlapping[nd]+1

if iw1[nd]<1
 isol1[nd]=isol1[nd]+iw1[nd]-1
 iw1[nd]=1
end

iw2[nd]=ij[nd]+stepsize[nd]+overlapping[nd]-1;

if iw2[nd]>size(mask)[nd]
iw2[nd]=size(mask)[nd]
end


isol2[nd]=isol1[nd]+stepsize[nd]-1

if isol2[nd]> iw2[nd]-iw1[nd]+1
 isol2[nd]=iw2[nd]-iw1[nd]+1
end

istore1[nd]=ij[nd];
istore2[nd]=istore1[nd]+isol2[nd]-isol1[nd];




end



warn("Test window $iw1 $iw2 $isol1 $isol2 $ij $istore1 $istore2 ")




#################################################
# Need to check how to work with aditional constraints...
#################################################



windowpoints=([iw1[i]:iw2[i] for i in 1:n]...);

# @time fw,s=divandrun(mask,pmn,xi,x,f,Labs,epsilon2; otherargs...)

fw,s=divandjog(mask[windowpoints...],([ x[windowpoints...] for x in pmn ]...),([ x[windowpoints...] for x in xi ]...),x,f,Labs,epsilon2; otherargs...)

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

# ndlast=sum(1-s.obsout)

# end
# end
# end
end


return fi,s #,pmnw,xw,maskw,ndlast






# return pmnw,n,Lscales,overlapping,stepsize
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
