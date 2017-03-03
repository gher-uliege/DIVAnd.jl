"""
Compute a variational analysis of arbitrarily located observations to calculate an estimate of the optimal value of epsilon2

 = divand_cv(mask,pmn,xi,x,f,len,epsilon2,nl,ne,...);

Perform an n-dimensional variational analysis of the observations `f` located at
the coordinates `x`. The output `factors` represent multipliction factors applied to epsilon2 which have been tested and the cvvalues the corresponding cross validation values.

The epsilon2 provided should be close the real one as the tests will be performed around


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

* `epsilon2`: signal-to-noise ratio of observations (if epsilon2 is a scalar).
    The larger this value is, the closer is the field `fi` to the
    observation. If epsilon2 is a scalar, then R is 1/epsilon2 I, where R is the observation error covariance matrix). If epsilon2 is a vector, then R is diag(epsilon2) or if epsilon2 is a matrix (a matrix-like project), then R is equal to epsilon2.

* `nl`: number of testing points around the current value of l. One means an addition point on both sides of the current L. Zero is allowed and means the parameter is not optimised.

* `ne`: number of testing points around the current value of epsilon2. Zero is allowed as for nl

# Optional input arguments specified as keyword arguments also as for divand


# Output:

* `bestfactorl`: best estimate of the multipliocation factor to apply to len

* `bestfactore`: best estimate of the multipliocation factor to apply to epsilon2

* `cvvales` : the cross validation values calculated

* `factors` : the tested multiplication factors

* `cvinter` : interpolated cv values for final optimisation

* `linter` : values of the factors at which the interpolation was done (in log scale)

* `epsilon2inter` : values of the factors at which the interpolation was done (in log scale)

"""


function divand_cv(mask,pmn,xi,x,f,len,epsilon2,nl,ne,method=0; otherargs...)

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
        logfactorsl=collect(linspace(-worderl,worderl,2*nl+1));
    else
        logfactorsl=[0]
    end
    factorsl=10.^logfactorsl;


    if ne>0
        logfactorse=collect(linspace(-wordere,wordere,2*ne+1));
    else
        logfactorse=[0]
    end
    factorse=10.^logfactorse;

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


            fi,s =  divandrun(mask,pmn,xi,x,f,len.*factorsl[i],epsilon2.*factorse[j]; otherargs...);
            residual=divand_residualobs(s,fi);
            nrealdata=sum(1-s.obsout);
            d0d=dot((1-s.obsout).*(s.yo),(s.yo));
            d0dmd1d=dot((1-s.obsout).*residual,(s.yo));

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
                cvval=divand_cvestimator(s,residual./(1-divand_diagHKobs(s)));
                epsilon2in[ip] = 1/5000;
            end



            if mymethod==2
                # alternate version to test: sampling
                # find(x -> x == 3,z)
                #   onsea=find(x->x == 0,s.obsout);
                onsea=find(s.obsout.==0);
                lonsea=length(onsea)
                #   warn("So",lonsea)
                # if optimisation is to be used, make sure to use the same reference random points
                srand(nrealdata)
                # otherwise you add noise to the cv field
                indexlist1=unique(collect(rand(1:lonsea,50*samplesforHK)))[1:samplesforHK]
                srand()
                indexlist=onsea[indexlist1];
                #   indexlist=collect(1:lonsea);
                residualc=zeros(length(residual));
                residualc[indexlist]=residual[indexlist]./(1-divand_diagHKobs(s,indexlist))
                scalefac=float(nrealdata)/float(samplesforHK)
                cvval=scalefac*divand_cvestimator(s,residualc)
                epsilon2in[ip] = 1/5000;
            end

            if mymethod==3
                work=(1-divand_GCVKiiobs(s));
                cvval=divand_cvestimator(s,residual./work);
                epsilon2in[ip] = 1/2000/work^2;
            end

            cvvalues[ip]=cvval;

        end
    end





    #####################


    # When no sampling is requested, just return CV value (for divand_qc) and a
    # multiplication factor for epsilon2 based on the Derozier adaptation idea
    # ll1= d0d/(d0d-d1d)-1
    #
    if (ne==0 && nl==0)
        warn("There is no parameter optimisation done ", nl, ne)
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





    # Now interpolate and find minimum using 1D divand or 2D divand depending on the situation




    if nl==0

        # interpolate only on epsilon
        epsilon2inter=collect(linspace(-wordere*1.1,1.1*wordere,101))
        maskcv = trues(size(epsilon2inter))
        pmcv = ones(size(epsilon2inter)) / (epsilon2inter[2]-epsilon2inter[1])
        lenin = wordere;

        m = Int(ceil(1+1/2))
        alpha = [binomial(m,k) for k = 0:m];
        alpha[1]=0;

        cvinter,scv = divandrun(maskcv,(pmcv,),(epsilon2inter,),(logfactorse,),cvvalues,lenin,epsilon2in;alpha=alpha,alphabc=0)


        bestvalue=findmin(cvinter)
        posbestfactor=bestvalue[2]
        cvval=bestvalue[1]
        bestfactor=10^epsilon2inter[posbestfactor]
        return bestfactor, cvval,cvvalues, logfactorse,cvinter,epsilon2inter

    end

if ne==0

    # interpolate only on L
    linter=collect(linspace(-worderl*1.1,1.1*worderl,101))
    maskcv = trues(size(linter))
    pmcv = ones(size(linter)) / (linter[2]-linter[1])
    lenin = worderl;

    m = Int(ceil(1+1/2))
    alpha = [binomial(m,k) for k = 0:m];
    alpha[1]=0;

    cvinter,scv = divandrun(maskcv,(pmcv,),(linter,),(logfactorsl,),cvvalues,lenin,epsilon2in;alpha=alpha,alphabc=0)


    bestvalue=findmin(cvinter)
    posbestfactor=bestvalue[2]
    cvval=bestvalue[1]
    bestfactor=10^linter[posbestfactor]
    return bestfactor, cvval,cvvalues, logfactorsl,cvinter,linter


end


# Otherwise 2D

xi2D,yi2D   = ndgrid(linspace(-worderl*1.1,1.1*worderl,71),linspace(-wordere*1.1,1.1*wordere,71))


maskcv = trues(xi2D)

# this problem has a simple cartesian metric
# pm is the inverse of the resolution along the 1st dimension
# pn is the inverse of the resolution along the 2nd dimension
pm2D = ones(xi2D) / (xi2D[2,1]-xi2D[1,1]);
pn2D = ones(xi2D) / (yi2D[1,2]-yi2D[1,1]);



# correlation length
lenin = (worderl,wordere);

m = Int(ceil(1+2/2))
alpha = [binomial(m,k) for k = 0:m];
alpha[1]=0;

cvinter,scv = divandrun(maskcv,(pm2D,pn2D),(xi2D,yi2D),(x2Ddata,y2Ddata),cvvalues,lenin,epsilon2in;alpha=alpha,alphabc=0)

bestvalue=findmin(cvinter)
posbestfactor=bestvalue[2]
cvval=bestvalue[1]
bestfactorl=10^xi2D[posbestfactor]
bestfactore=10^yi2D[posbestfactor]
return bestfactorl,bestfactore, cvval,cvvalues, x2Ddata,y2Ddata,cvinter,xi2D,yi2D






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


