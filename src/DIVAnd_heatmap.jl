"""
Computes a  heatmap based on locations of observations using kernel density estimation (probability density field whose integral over the domain is one)

dens,Ltuple,LCV,LSCV = DIVAnd_heatmap(mask,pmn,xi,x,inflation,Labs;Ladaptiveiterations=0,myheatmapmethod="DataKernel",
    optimizeheat=true,nmax=1000,otherargs...)


# Input:
*  `mask`: mask as usual
*  `pmn` : tuple of metrics as usual
*  `xi`: tuple of coordinates of the grid for the heatmap
*  `x` : tuple of coordinates of observations
*  `inflation`: array generally of ones. For some applications an observation can carry a different weight which is then encoded in the array
*  `Labs` : the length scales for DIVAnd. Here their meaning is the spread (bandwidth) of the observations for the Kernel calculation
*              if zero is provided, the routine applies an empirical estimate, returned in the Ltuple output.

*   `Ladaptiveiterations`: adaptive scaling where the length scales are adapted on the data density already estimated. You can iterate. Default "0"
*   `optimizeheat` : boolean which can turn on or off an algorithmic optimisation. Results should be identical. Default is to optimize
*   `myheatmapmethod`: can be "Automatic", "GridKernel" or "DataKernel" (Results should be very similar except near boundaries)
*   `nmax`: maximum number of data points. If actual data size is larger, approximatively nmax superobservations are calculated and a warning issued.

*   `otherargs...`: all other optional arguments DIVAndrun can take (advection etc)

# Output:
*  `dens`: data density field (integral is one)
*  `Ltuple` : The bandwidthth used (either the input value or the calculated ones)
*  `LCV` : Likelihood Cross Validation estimator value (the higher the better) leave one out approach
*  `LSCV` : Least Square Cross Validation estimator (the lower the better) leave one out approach
"""
function DIVAnd_heatmap(
    mask,
    pmn,
    xi,
    xin,
    inflationin,
    Labs;
    Ladaptiveiterations = 0,
    myheatmapmethod = "DataKernel",
    optimizeheat = true,
    nmax = 10000,
    otherargs...,
)
#
    x = xin
    inflation = inflationin
    idx = []
    if size(inflationin)[1] > nmax
        @warn "Data array size $(size(inflationin)) is larger then maximum $nmax. Superobservations will be created. To avoid, increase nmax to the desired number of superobs"

        x, inflation, sumw, varp, idx = DIVAnd_superobs(
            xin,
            ones(size(inflationin)),
            nmax;
            weights = inflationin,
            intensive = false,
        )
       #@show size(inflationin),nmax,size(inflation)
    end



# Create output array on the same grid as mask pmn and xi
    dens2 = zeros(Float64, size(mask))
# Dimensionality of the problem and number of points
    DIMS = ndims(mask)
    NP = size(inflation)[1]
    NPI = sum(inflation)
    selfvalue = zeros(Float64, NP)
    svf = statevector_init((mask,))
    LCV = 0.0
    LSCV = 0.0

    LHEAT = Labs

    mymethod = myheatmapmethod

    trytooptimize = optimizeheat

# If automatic selection take the one with the lower number of covariances to calculate
    if myheatmapmethod == "Automatic"

        mymethod = "DataKernel"

        if NP > sum(mask .== true)
            mymethod = "GridKernel"
        end


    end

    if mymethod == "GridKernel"

        @warn "Method GridKernel does not allow for cross validation. If you need the latter, force myheatmapmethod to DataKernel"

    end

    if trytooptimize == false
        @warn "Unoptimized versions do not allow for cross validation"
    end


#
    if Labs == 0
        # Estimate
          # Empirial estimate Silverman's (1986) rule of thumb

        varx = zeros(Float64, DIMS)
          #varxb=zeros(Float64,DIMS)
        LF = zeros(Float64, DIMS)
        #fromgausstodivaL=[0.8099,0.62,0.62,0.62,0.62,0.62,0.62]
        for i = 1:DIMS
            meanxo = sum(inflation .* x[i]) / sum(inflation)
            varx[i] = sum(inflation .* (x[i] .- meanxo).^2) / sum(inflation)
               #meanxob=sum(x[i])/size(inflation)[1]
            #varxb[i]=sum((x[i].-meanxob).^2)/size(inflation)[1]
               #@show i,meanxo,varx[i],meanxob,varxb[i]
            LF[i] = sqrt(varx[i]) / ((DIMS + 2.0) * NPI / 4.0)^(1.0 / (4.0 + DIMS))
               #@show LF[i]
               #LF[i]=sqrt(varxb[i])/((DIMS+2.0)*NP/4.0)^(1.0/(4.0+DIMS))
               #@show LF[i]

        end

        if idx != []
          #@show "was binned"
          #@show LF.*idx
            for i = 1:DIMS
                if LF[i] * idx[i] < 1
                    @warn "Superobserving made length scale $(LF[i]) in direction $(i) estimate too small, forced to bin size $(1/idx[i])"
                    LF[i] = 1.0 / idx[i]
                end
            end
        end

        LHEAT = (LF...,)
        #@show "Estimation", LHEAT

    end



    #
    Lfortuple = Array{Any}(undef, DIMS)
    Ltuple = LHEAT


#
    inflationsum = 0





    for Literations = 1:1+Ladaptiveiterations

        inflationsum = 0
        dens2 .= 0.0


        xxx = Array{Any}(undef, DIMS)
        if trytooptimize
        #@show "Try to calculate a decomposition"
            #Decompose once and for all
               #if mymethod=="DataKernel"
            #FIopt,Sopt=DIVAnd.DIVAndrun(mask,pmn,xi,x,inflation,Ltuple, 1.0E10 ;otherargs...)
            #svf=statevector_init((mask,))

               #end
               #if mymethod=="GridKernel"
            FIopt, Sopt = DIVAnd.DIVAndrun(
                mask,
                pmn,
                xi,
                x,
                inflation,
                Ltuple,
                1.0E10;
                otherargs...,
            )

               #end
        end

# VERSION A: covariance of one data point with grid points
        if mymethod == "DataKernel"

            for myi = 1:NP

                if trytooptimize
                    eiarr = zeros(size(inflation))
                    eiarr[myi] = 1
                    work1 = Sopt.H' * eiarr
                    vv = Sopt.P.factors.PtL \ work1
                    vb = Sopt.P.factors.UP \ vv
                   # Possible place for slight performance improvement. Do the integral in state-space with volumes created before.
                    FI, = statevector_unpack(svf, vb)
                    integ = DIVAnd_integral(mask, pmn, FI)
                    selfvalue[myi] = (work1 ⋅ vb) / integ

                else


                    for ii = 1:DIMS
                        xxx[ii] = [x[ii][myi]]
                    end
            #@show xxx,LF,size(xi[1]),size(pmn[1])
            #  Use of WOODBURY and decomposed B in s could make it faster, done in the optimized version
                    FI, S = DIVAnd.DIVAndrun(
                        mask,
                        pmn,
                        xi,
                        (xxx...,),
                        [1.0],
                        Ltuple,
                        0.001;
                        otherargs...,
                    )

    # Add here the constraint that each integral is one: accepts non unit values vi inflation to reflect mulitple observations
    # Also accept errors on obs? But how ? In the integral so it has less influence on overall sum (which will be scaled again?)
    #
                    integ = DIVAnd_integral(mask, pmn, FI)
        #@show integ

                end
                if integ != 0
                    dens2 = dens2 .+ inflation[myi] * FI / integ
                    inflationsum = inflationsum + inflation[myi]
            #else
            #    @show "?? Not on grid ?",integ,sum(FI),xxx
                end



            end


            dens2 = dens2 / inflationsum

            if trytooptimize
                dens2x = statevector_pack(svf, (dens2,))
            #@show sum(dens2),DIVAnd_integral(mask,pmn,dens2),NP,inflationsum
                selfvalueerr = (inflationsum .* Sopt.H * dens2x .-
                                inflation .* selfvalue) ./ (inflationsum .- inflation)
                selfvalueerr[selfvalueerr.<0.0] .= 0.0

                logvalueerr = log.(selfvalueerr)


                finalweights = inflation
               #or
               #finalweights=ones(Float64,size(inflation))

                finalweights[isnan.(selfvalueerr)] .= 0
                finalsum = sum(finalweights)
                selfvalueerr[isnan.(selfvalueerr)] .= 0
                logvalueerr[isnan.(logvalueerr)] .= 0
            #errestim=inflation.*(selfvalue-Sopt.H*dens2x)./(inflationsum.-inflation)
               #@show DIVAnd_integral(mask,pmn,dens2.^2), DIVAnd_integral(mask,pmn,dens2.^2)-2*sum(selfvalueerr)/inflationsum,sum(logvalueerr)

                LCV = sum(finalweights .* logvalueerr) / finalsum
                LSCV = DIVAnd_integral(mask, pmn, dens2.^2) -
                       2 * sum(finalweights .* selfvalueerr) / finalsum
               #or

               #LCV=sum(inflation.*logvalueerr)/inflationsum
               #LSCV=DIVAnd_integral(mask,pmn,dens2.^2)-2*sum(inflation.*selfvalueerr)/inflationsum

            end
        end


# VERSION B: covariance of one grid point with all data points
        if mymethod == "GridKernel"

            if trytooptimize
                svsize = sum(mask .== true)
                xdens = zeros(Float64, svsize)
                xval = zeros(Float64, svsize)
              #@show "OPti",svsize
                for myi = 1:svsize
                    eiarr = zeros(Float64, svsize)
                    eiarr[myi] = 1.0
                    vv = Sopt.P.factors.PtL \ eiarr
                    vb = Sopt.P.factors.UP \ vv
                    FI, = statevector_unpack(svf, vb)
                    integ = DIVAnd_integral(mask, pmn, FI)
                    vb = vb / integ
                    xdens[myi] = sum(inflation .* (Sopt.H * vb))
                end
                dens2, = statevector_unpack(svf, xdens)
                dens2 = dens2 / DIVAnd_integral(mask, pmn, dens2)

            else

                xaugmented = Array{Any}(undef, DIMS)
                Rinf = deepcopy(inflation)
                Rinf .= 1.0E10
                Raugmented = [Rinf..., 0.000001]
                valaugmented = [inflation..., 1.0]
                for myi in eachindex(mask)

                    if mask[myi]
                    # Only on grid:
                    # Add one virtual data point which is the grid point.
                        for ii = 1:DIMS
                            xaugmented[ii] = [x[ii]..., xi[ii][myi]]
                        end

                    # The real data points with infinite R
                    # Use of WOODBURY and decomposed B in s could make it faster
                        FI, s = DIVAnd.DIVAndrun(
                            mask,
                            pmn,
                            xi,
                            (xaugmented...,),
                            valaugmented,
                            Ltuple,
                            Raugmented;
                            otherargs...,
                        )

                    # Calculate normalization constant
                        integ = DIVAnd_integral(mask, pmn, FI)
                    # After analysis, retrieve via S the analysis at those points and apply normalization constant
                    # as well as inflation
                        bidon = (s.obsconstrain.H) * statevector_pack(s.sv, (FI,))
                        dens2[myi] = sum(bidon[1:end-1] .* inflation) / integ

                    #@show myi,integ


                    end

                end

                dens2 = dens2 / DIVAnd_integral(mask, pmn, dens2)
            end
        end








        # Optional L scaling

        if Ladaptiveiterations > 0
            #@show Literations,size(Ltuple[1])
            lambda = DIVAnd_scaleL(mask, pmn, dens2)
            for jj = 1:DIMS
                #@show LHEAT[jj]
                Lfortuple[jj] = LHEAT[jj] .* lambda
            end
            Ltuple = (Lfortuple...,)
        end




    end

    dens2[.!mask] .= NaN

    return dens2, Ltuple, LCV, LSCV

end
#

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


