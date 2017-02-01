"""
Compute a variational analysis of arbitrarily located observations to calculate the clever poor man's error

cpme = divand_aexerr(mask,pmn,xi,x,f,len,lambda,...);

Perform an n-dimensional variational analysis of the observations `f` located at
the coordinates `x`. The array `cpme` represent the error field at the grid
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

* `lambda`: signal-to-noise ratio of observations (if lambda is a scalar).
    The larger this value is, the closer is the field `fi` to the
    observation. If lambda is a scalar, then R is 1/lambda I, where R is the observation error covariance matrix). If lambda is a vector, then R is diag(lambda) or if lambda is a matrix (a matrix-like project), then R is equal to lambda.


# Optional input arguments specified as keyword arguments also as for divand


# Output:

* `aexerr`: the almost exact error
"""


function divand_aexerr(mask,pmn,xi,x,f,len,lambda; otherargs...)

# Hardwired value:

finesse=2

# No need to make an approximation if it is close to cost of direct calculation
upperlimit=0.4

# check inputs

if !any(mask[:])
  error("no sea points in mask");
end

# Decide how many additional fake points are needed with almost zero weight


# Assuming uniform grid
n = ndims(mask)
nsamp=ones(n);
npgrid=1;
npneeded=1;

for i=1:n
	if isa(len,Number)
		Labs = len;
	elseif isa(len,Tuple)

		if isa(len[1],Number)
			Labs = len[i]
			else
			Labs=len[i][1]
		end

	end
	npgrid=npgrid*size(mask)[i];
	nsamp[i]=Labs*pmn[i][1]/finesse;
	npneeded=npneeded*size(mask)[i]/nsamp[i];
end


ndata=size(f)[1];


if npneeded>upperlimit*npgrid
# 
   return 0,0,0,0

end


npongrid=Int(ceil(maximum([npgrid/10^n,npneeded-ndata])));

randindexes=ones(Int,n,npongrid);

for i=1:n
randindexes[i,:]=rand(1:size(mask)[i],npongrid);
end


# add npongrind fake points onto the grid with zero value and very high R value
xfake=x;
#for i=1:n
#xfake[i]=append!(x[i], xi[i][randindexes[i,:]])
#end

# Make an analysis with those fake points and very low snr to get B at those locations 

xfake=x;
ffake=f;
lambdafake=0.0001;
f1,s1=divandrun(mask,pmn,xi,xfake,ffake,len,lambdafake; otherargs...);

# Interpolate B on the final grid with high snr
Batdatapoints=divand_erroratdatapoints(s1);

# Would be nice to use semi norm here ...
m = Int(ceil(1+n/2))
  # alpha is the (m+1)th row of the Pascal triangle:
  # m=0         1
  # m=1       1   1
  # m=1     1   2   1
  # m=2   1   3   3   1
  # ...
alpha = [binomial(m,k) for k = 0:m];
alpha[1]=0;

Bjmb,s1=divandrun(mask,pmn,xi,xfake,Batdatapoints,len*2,10)#; alpha=alpha )#), otherargs...)

# Now do the same with normal snr to get real error at the "data" points
lambdafake=lambda
f1,s1=divandrun(mask,pmn,xi,xfake,ffake,len,lambdafake; otherargs...);
Errdatapoints=divand_erroratdatapoints(s1);


# Now get error reduction terms
ffake=Batdatapoints-Errdatapoints;

# Interpolate error reduction term 
f1,s1=divandrun(mask,pmn,xi,xfake,ffake,len./1.70766,100; otherargs...);

# Calculate final error
aexerr=Bjmb-f1;





# Provide the error field, the background field for additional scaling and the analysis itself

return npongrid,randindexes,aexerr,Bjmb,x,xi,randindexes

#return aexerr,berr,fi,s

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


