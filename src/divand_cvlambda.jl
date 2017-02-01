"""
Compute a variational analysis of arbitrarily located observations to calculate an estimate of the optimal value of lambda

bestfactor,cvvalues,factors, cvinter,laminter = divand_cvlambda(mask,pmn,xi,x,f,len,lambda,...);

Perform an n-dimensional variational analysis of the observations `f` located at
the coordinates `x`. The output `factors` represent multipliction factors applied to lambda which have been tested and the cvvalues the corresponding cross validation values.

The lambda provided should be close the real one as the tests will be performed around


The analysus is defined by the coordinates `xi` and the scales factors `pmn`.

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

* `lambda`: signal-to-noise ratio of observations (if lambda is a scalar).
    The larger this value is, the closer is the field `fi` to the
    observation. If lambda is a scalar, then R is 1/lambda I, where R is the observation error covariance matrix). If lambda is a vector, then R is diag(lambda) or if lambda is a matrix (a matrix-like project), then R is equal to lambda.


# Optional input arguments specified as keyword arguments also as for divand


# Output:

* `bestfactor`: best estimate of the multipliocation factor to apply to lambda

* `cvvales` : the cross validation values calculated 

* `factors` : the tested multiplication factors

* `cvinter` : interpolated cv values for final optimisation 

* `laminter` : values of the factors at which the interpolation was done (in log scale)

"""


function divand_cvlambda(mask,pmn,xi,x,f,len,lambda; otherargs...)

# check inputs

if !any(mask[:])
  error("no sea points in mask");
end
# For the moment, hardwired values
switchvalue=100;
worder=1.5;
nsamp=2;

# sample multiplication factor to optimise in log space
logfactors=collect(linspace(-worder,worder,2*nsamp+1));
factors=10.^logfactors;
cvvalues=0.*factors;




for i=1:size(factors)[1]

fi,s =  divandrun(mask,pmn,xi,x,f,len,lambda.*factors[i]; otherargs...);
residual=divand_residualobs(s,fi);
nrealdata=sum(1-s.obsout);


if nrealdata<switchvalue
   cvval=divand_cvestimator(s,residual./(1-divand_diagHKobs(s)));
         else
   cvval=divand_cvestimator(s,residual./(1-divand_GCVKiiobs(s)));	 
end

cvvalues[i]=cvval;

end
# Now interpolate and find minimum using 1D divand

laminter=collect(linspace(-worder*1.1,1.1*worder,101))

maskcv = trues(size(laminter))

# this problem has a simple cartesian metric
# pm is the inverse of the resolution along the 1st dimension
# pn is the inverse of the resolution along the 2nd dimension

pmcv = ones(size(laminter)) / (laminter[2]-laminter[1])


# correlation length
lenin = worder;

# signal-to-noise ratio
lambdain = 500;

# fi is the interpolated field
# TODO adapt for seminorm
m = Int(ceil(1+1/2))
  # alpha is the (m+1)th row of the Pascal triangle:
  # m=0         1
  # m=1       1   1
  # m=1     1   2   1
  # m=2   1   3   3   1
  # ...
alpha = [binomial(m,k) for k = 0:m];
alpha[1]=0;
cvinter,scv = divandrun(maskcv,(pmcv,),(laminter,),(logfactors,),cvvalues,lenin,lambdain;alpha=alpha)

posbestfactor=findmin(cvinter)[2]
bestfactor=10^laminter[posbestfactor]
return bestfactor, cvvalues, factors,cvinter,laminter

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


