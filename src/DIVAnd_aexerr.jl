"""
    aexerr,Bref,fa,sa = DIVAnd_aexerr(mask,pmn,xi,x,f,len,epsilon2;...);



# Input: same as for DIVAndrun

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

# Optional input arguments specified as keyword arguments also as for DIVAnd


# Output:

* `aexerr`: the almost exact error

* `Bref`: the background error for error scaling by background `aexerr./Bref`

* `fa`: the analysis (with low impact fake data): DO NOT USE UNLESS YOU KNOW WHAT YOU ARE DOING

* `sa`: the associated structure

Compute a variational analysis of arbitrarily located observations to calculate the almost exact error

"""
function DIVAnd_aexerr(mask,pmn,xi,x,f,len,epsilon2; otherargs...)



    # Hardwired value:
    finesse=1.7

    # No need to make an approximation if it is close to cost of direct calculation
    upperlimit=0.4




    oriR=DIVAnd_obscovar(epsilon2,size(f)[1]);

    epsilonref=mean(diag(oriR));

    epsilonslarge=maximum([1E6,epsilonref*1E6]);



    # Decide how many additional fake points are needed with almost zero weight


    #
    n = ndims(mask)
    nsamp=ones(n);
    npgrid=1;
    npneeded=1;
    Labspmnmin=zeros(n)

    for i=1:n
        if isa(len,Number)
            Labspmnmin[i] = len*minimum(pmn[i]);
        elseif isa(len,Tuple)

            if isa(len[1],Number)
                Labspmnmin[i] = len[i]*minimum(pmn[i]);

            else
                Labspmnmin[i] = minimum(len[i].*pmn[i])

            end

        end
        npgrid=npgrid*size(mask)[i];
        nsamp[i]=Labspmnmin[i]/finesse;
        npneeded=npneeded*size(mask)[i]/nsamp[i];

    end


    ndata=size(f)[1];


    if npneeded>upperlimit*npgrid
        # No need to make an approximation if it is close to cost of direct calculation
        return 0,0,0,0
        # Need to catch this event outside and use direct error calculation instead
    end


    npongrid=Int(ceil(maximum([npgrid/10^n,npneeded-ndata])));

    randindexes=ones(Int,npongrid);

    nsa=Int(ceil(npgrid/npongrid));

    randindexes=collect(1:nsa:npgrid);

    ncv=size(randindexes)[1];

    # add npongrind fake points onto the grid with zero value and very high R value

    ffake=deepcopy(f);

    ffake=append!(ffake, 0. * xi[1][randindexes]);
    Rfake=blkdiag(oriR,DIVAnd_obscovar(epsilonslarge,ncv));
    xfake=tuple([append!(copy(x[i]), xi[i][randindexes]) for i=1:n]...)

    # Make an analysis with those fake points and very low snr to get B at those locations


    epsilon2fake=10_000.
    f1,s1=DIVAndrun(mask,pmn,xi,xfake,ffake,len,epsilon2fake; otherargs...);


    # Interpolate B on the final grid with high snr

    # First get B, the error of the previous analysis with bad data at the data locations
    Batdatapoints=DIVAnd_erroratdatapoints(s1);

    # Now use semi norm here ...
    m = Int(ceil(1+n/2))
    alpha = [binomial(m,k) for k = 0:m];
    alpha[1]=0;


    # Analyse with semi-norm and larger length scales
    Bjmb,s1=DIVAndrun(mask,pmn,xi,xfake,Batdatapoints,len*2,1/20; alpha=alpha, otherargs...)

    Bjmb=max.(Bjmb,0)

    # Now do the same with normal snr to get real error at the "data" points
    # incidentally fa and sa are almost the real analysis
    # @show typeof(Rfake)
    # @show issubtype(typeof(Rfake),Union{AbstractArray{Float64,1},AbstractArray{Float64,2}})

    fa,sa=DIVAndrun(mask,pmn,xi,xfake,ffake,len,Rfake; otherargs...);
    Errdatapoints=DIVAnd_erroratdatapoints(sa);


    # Now get error reduction terms
    ffake=Batdatapoints-Errdatapoints;

    # Interpolate error reduction term
    # The factor 1.70677 is the best one in 2D but should be slightly different for other dimensions
    # Could be a small improvement. Also used in DIVAnd_cpme

    f1,s1=DIVAndrun(mask,pmn,xi,xfake,ffake,len./1.70766,1.0/100.0; otherargs...);

    # Calculate final error
    aexerr=max.(Bjmb-f1,0);





    # Provides the almost error field aexerr, the background field Bjmb for additional scaling and the analysis fa itself with its strucure sa

    return aexerr,Bjmb,fa,sa



end

# Copyright (C) 2008-2019 Alexander Barth <barth.alexander@gmail.com>
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

# LocalWords:  fi DIVAnd pmn len diag CovarParam vel ceil moddim fracdim
