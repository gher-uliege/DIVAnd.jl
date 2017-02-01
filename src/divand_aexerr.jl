"""
Compute a variational analysis of arbitrarily located observations to calculate the clever poor man's error

cpme = divand_aexerr(mask,pmn,xi,x,f,len,epsilon2,...);

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

* `epsilon2`: error variance of the observations (normalized by the error variance of the background field). `epsilon2` can be a scalar (all observations have the same error variance and their errors are decorrelated), a vector (all observations can have a difference error variance and their errors are decorrelated) or a matrix (all observations can have a difference error variance and their errors can be correlated). If `epsilon2` is a scalar, it is thus the *inverse of the signal-to-noise ratio*.

# Optional input arguments specified as keyword arguments also as for divand


# Output:

* `aexerr`: the almost exact error

* `Bjmb`: the background error

* `fa`: the analysis (with low impact fake data)

* `sa`: the associated structure
"""


function divand_aexerr(mask,pmn,xi,x,f,len,epsilon2; otherargs...)



# Hardwired value:

finesse=1.7

# No need to make an approximation if it is close to cost of direct calculation
upperlimit=0.4

# check inputs

if !any(mask[:])
  error("no sea points in mask");
end


oriR=divand_obscovar(epsilon2,size(f)[1]);

epsilonref=mean(diag(oriR));

epsilonslarge=maximum([1E6,epsilonref*1E6]);



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

randindexes=ones(Int,npongrid);

nsa=Int(ceil(npgrid/npongrid));

#randindexes=rand(1:npgrid,npongrid);
randindexes=collect(1:nsa:npgrid);

ncv=size(randindexes)[1];

# add npongrind fake points onto the grid with zero value and very high R value
#xfake=x;
#for i=1:n
# xfake[i]=append!(x[i], xi[i][randindexes])
#end
ffake=deepcopy(f);
ffake=append!(ffake, 0.*xi[1][randindexes]);

Rfake=blkdiag(oriR,divand_obscovar(epsilonslarge,ncv));

#xcc=deepcopy(x);
xfake=tuple([append!(copy(x[i]), xi[i][randindexes]) for i=1:n]...)

# Make an analysis with those fake points and very low snr to get B at those locations 
#xfake=x;
#ffake=f;

epsilon2fake=10_000;
f1,s1=divandrun(mask,pmn,xi,xfake,ffake,len,epsilon2fake; otherargs...);

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

Bjmb,s1=divandrun(mask,pmn,xi,xfake,Batdatapoints,len*2,1/20; alpha=alpha, otherargs...)



# Now do the same with normal snr to get real error at the "data" points
# incidentally fa and sa are almost the real analysis
fa,sa=divandrun(mask,pmn,xi,xfake,ffake,len,Rfake; otherargs...);
Errdatapoints=divand_erroratdatapoints(sa);


# Now get error reduction terms
ffake=Batdatapoints-Errdatapoints;

# Interpolate error reduction term 
f1,s1=divandrun(mask,pmn,xi,xfake,ffake,len./1.70766,1/100; otherargs...);

# Calculate final error
aexerr=Bjmb-f1;





# Provide the error field, the background field for additional scaling and the analysis itself

return aexerr,Bjmb,fa,sa

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


