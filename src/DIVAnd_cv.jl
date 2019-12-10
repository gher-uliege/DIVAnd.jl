# nofmt

# disable formatting: https://github.com/domluna/JuliaFormatter.jl/issues/72


"""

    bestfactorl,bestfactore, cvval,cvvalues, x2Ddata,y2Ddata,cvinter,xi2D,yi2D = DIVAnd_cv(mask,pmn,xi,x,f,len,epsilon2,nl,ne,method;...);

Performs a cross validation to estimate the analysis parameters
(correlation length and signal-to-noise ratio).

# Input

Same as for `DIVAndrun` with three more parameters `nl`,`ne` and `method`

* `mask`: binary mask delimiting the domain. true is inside and false outside.
For oceanographic application, this is the land-sea mask.

* `pmn`: scale factor of the grid. pmn is a tuple with n elements. Every
       element represents the scale factor of the corresponding dimension. Its
       inverse is the local resolution of the grid in a particular dimension.

* `xi`: tuple with n elements. Every element represents a coordinate
  of the final grid on which the observations are interpolated

* `x`: tuple with n elements. Every element represents a coordinate of
  the observations

* `f`: value of the observations *minus* the background estimate (m-by-1 array).
    (see note)

* `len`: correlation length

* `epsilon2`: error variance of the observations (normalized by the error variance of the background field). `epsilon2` can be a scalar (all observations have the same error variance and their errors are decorrelated), a vector (all observations can have a difference error variance and their errors are decorrelated) or a matrix (all observations can have a difference error variance and their errors can be correlated). If `epsilon2` is a scalar, it is thus the *inverse of the signal-to-noise ratio*.

* `nl`: number of testing points around the current value of L. `1` means one additional point on both sides of the current L. `0` is allowed and means the parameter is not optimised.

* `ne`: number of testing points around the current value of epsilon2. `0` is allowed as for `nl`

* `method`: cross validation estimator method
  1: full CV
  2: sampled CV
  3: GCV
  0: automatic choice between the three possible ones, default value

* Optional input arguments specified via keyword arguments are the same as for `DIVAnd`


# Output:

* `bestfactorl`: best estimate of the multiplication factor to apply to len

* `bestfactore`: best estimate of the multiplication factor to apply to epsilon2

* `cvvales` : the cross validation values calculated

* `factors` : the tested multiplication factors

* `cvinter` : interpolated cv values for final optimisation

* `X2Data, Y2Data` : coordinates of sampled cross validation in `L,epsilon2` space . Normally only used for debugging or plotting

* `Xi2D, Yi2D` : coordinates of interpolated estimator . Normally only used for debugging or plotting




The output `bestfactorl` and `bestfactore` represent multiplication factors which should be applied to `L` and `epsilon2`.


The `len` and `epsilon2` provided should be close the real one as the tests will be performed around.


"""
function DIVAnd_cv(mask,pmn,xi,x,f,len,epsilon2,nl,ne,method=0; otherargs...)

    # check inputs


    # For the moment, hardwired values
    switchvalue1=130;
    switchvalue2=1000;
    samplesforHK=75;
    # with of window so sample in log scale, so order of magnitudes worder
    # For length scales, one order of magnitude; make sure the grid is fine enough
    worderl=1.;
    # For noise, almost two order of magnitutes
    wordere=1.5;


    # sample multiplication factor to optimise in log space
    if nl>0
        logfactorsl=collect(range(-worderl,stop=worderl,length=2*nl+1));
    else
        logfactorsl=[0]
    end
    factorsl=10 .^ logfactorsl;


    if ne>0
        logfactorse=collect(range(-wordere,stop=wordere,length=2*ne+1));
    else
        logfactorse=[0]
    end
    factorse=10 .^ logfactorse;

    # cvvalues at the locations
    cvvalues=zeros((2*nl+1)*(2*ne+1));

    # For later interpolation quality might vary
    epsilon2in=zeros((2*nl+1)*(2*ne+1));

    x2Ddata=zeros((2*nl+1)*(2*ne+1));
    y2Ddata=zeros((2*nl+1)*(2*ne+1));

    # Define method used
    # 1: full CV
    # 2: sampled CV
    # 3: GCV
    # 0: automatic choice, default value
    # For automatic choice this will be done later once the exact number of useful data is known


    d0d=0
    d0dmd1d=0
    ip=0
    for i=1:size(factorsl)[1]
        for j=1:size(factorse)[1]

            ip=ip+1
            x2Ddata[ip]=logfactorsl[i];
            y2Ddata[ip]=logfactorse[j];


            fi,s =  DIVAndrun(mask,pmn,xi,x,f,len.*factorsl[i],epsilon2.*factorse[j]; otherargs...);
            residual = DIVAnd_residualobs(s,fi);
            obsin = .!s.obsout
            nrealdata = sum(obsin)
            d0d = s.obsconstrain.yo[obsin] ⋅ s.obsconstrain.yo[obsin]
            d0dmd1d = s.obsconstrain.yo[obsin] ⋅ residual[obsin]

            # Determine which method to use

            if method==0

                mymethod=2
                if nrealdata < switchvalue1
                    mymethod=1
                end
                if nrealdata > switchvalue2
                    mymethod=3
                end


            else
                mymethod=method
            end

            # Check if more samples than data are asked switch to direct method

            if mymethod==2
                if nrealdata < samplesforHK
                    mymethod=1
                end
            end




            # TO DO : THINK ABOUT A VERSION WITH 30 REAL ESTIMATES OF KII AND THE RESIDUAL ONLY THERE
            # c
            # unique(collect(rand(1:1000,200)))[1:30]

            if mymethod==1
                cvval = DIVAnd_cvestimator(s,residual ./ (1 .- DIVAnd_diagHKobs(s)));
                epsilon2in[ip] = 1/5000;
            end



            if mymethod==2
                # alternate version to test: sampling
                # find(x -> x == 3,z)
                #   onsea=find(x->x == 0,s.obsout);
                onsea = findall(s.obsout.==0);
                lonsea=length(onsea)
                #   @warn "So",lonsea
                # if optimisation is to be used, make sure to use the same reference random points
                Random.seed!(nrealdata)

                # otherwise you add noise to the cv field
                indexlist1=unique(collect(rand(1:lonsea,50*samplesforHK)))[1:samplesforHK]
                Random.seed!()

                indexlist=onsea[indexlist1];
                #   indexlist=collect(1:lonsea);
                residualc=zeros(length(residual));
                residualc[indexlist]=residual[indexlist]./(1 .- DIVAnd_diagHKobs(s,indexlist))
                scalefac=float(nrealdata)/float(samplesforHK)
                cvval=scalefac*DIVAnd_cvestimator(s,residualc)
                epsilon2in[ip] = 1/5000;
            end

            if mymethod==3
                work=(1-DIVAnd_GCVKiiobs(s));
                cvval=DIVAnd_cvestimator(s,residual./work);
                epsilon2in[ip] = 1/2000/work^2;
            end

            cvvalues[ip]=cvval;

        end
    end





    #####################


    # When no sampling is requested, just return CV value (for DIVAnd_qc) and a
    # multiplication factor for epsilon2 based on the Derozier adaptation idea
    # ll1= d0d/(d0d-d1d)-1
    #
    if (ne==0 && nl==0)
        @warn "There is no parameter optimisation done (nl=$nl, ne=$ne)"
        ll1= d0d/(d0dmd1d)-1;
        eps1=1/ll1;
        if ndims(epsilon2) == 0
            eps2=epsilon2;
        elseif ndims(epsilon2) == 1
            eps2 = mean(epsilon2);
        else
            eps2 = mean(diag(epsilon2));
        end
        return cvvalues[1],eps1/eps2
    end





    # Now interpolate and find minimum using 1D DIVAnd or 2D DIVAnd depending on the situation




    if nl==0

        # interpolate only on epsilon
        epsilon2inter=collect(range(-wordere*1.1,stop=1.1*wordere,length=101))
        maskcv = trues(size(epsilon2inter))
        pmcv = ones(size(epsilon2inter)) / (epsilon2inter[2]-epsilon2inter[1])
        lenin = wordere;

        m = Int(ceil(1+1/2))
        alpha = [binomial(m,k) for k = 0:m];
        alpha[1]=0;

        cvinter,scv = DIVAndrun(maskcv,(pmcv,),(epsilon2inter,),(logfactorse,),cvvalues,lenin,epsilon2in;alpha=alpha,alphabc=0)


        bestvalue=findmin(cvinter)
        posbestfactor=bestvalue[2]
        cvval=bestvalue[1]
        bestfactor=10^epsilon2inter[posbestfactor]
        return bestfactor, cvval,cvvalues, logfactorse,cvinter,epsilon2inter

    end

if ne==0

    # interpolate only on L
    linter=collect(range(-worderl*1.1,stop=1.1*worderl,length=101))
    maskcv = trues(size(linter))
    pmcv = ones(size(linter)) / (linter[2]-linter[1])
    lenin = worderl;

    m = Int(ceil(1+1/2))
    alpha = [binomial(m,k) for k = 0:m];
    alpha[1]=0;

    cvinter,scv = DIVAndrun(maskcv,(pmcv,),(linter,),(logfactorsl,),cvvalues,lenin,epsilon2in;alpha=alpha,alphabc=0)


    bestvalue=findmin(cvinter)
    posbestfactor=bestvalue[2]
    cvval=bestvalue[1]
    bestfactor=10^linter[posbestfactor]
    return bestfactor, cvval,cvvalues, logfactorsl,cvinter,linter


end


# Otherwise 2D

maskcv,(pm2D,pn2D),(xi2D,yi2D) = DIVAnd_rectdom(
    range(-worderl*1.1,stop=worderl*1.1,length=71),
    range(-wordere*1.1,stop=wordere*1.1,length=71))

# correlation length
lenin = (worderl,wordere);

m = Int(ceil(1+2/2))
alpha = [binomial(m,k) for k = 0:m];
alpha[1]=0;

cvinter,scv = DIVAndrun(maskcv,(pm2D,pn2D),(xi2D,yi2D),(x2Ddata,y2Ddata),cvvalues,lenin,epsilon2in;alpha=alpha,alphabc=0)

bestvalue=findmin(cvinter)
posbestfactor=bestvalue[2]
cvval=bestvalue[1]
bestfactorl=10^xi2D[posbestfactor]
bestfactore=10^yi2D[posbestfactor]
return bestfactorl,bestfactore, cvval,cvvalues, x2Ddata,y2Ddata,cvinter,xi2D,yi2D






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
