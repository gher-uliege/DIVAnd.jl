"""


stepsize,overlapping,isdirect = DIVAnd_fittocpu(Lpmnrange,gridsize,latercsteps,moddim=[]);

# Creates a list of windows for subsequent domain decomposition
# Also calculates already the subsampling steps csteps for the preconditionners

# Input:

* `Lpmnrange`: For every dimension the minimum and maximum correlation length scaled by the local resolution (i.e. the product between L and pm (pn,...))
* `gridsize`: number of points in each direction (size(mask))
* `latercsteps`:  coarsening steps used later if a lower resolution model is used for preconditioning. 
* `moddim`: modulo for cyclic dimension (vector with n elements). Zero is used for non-cyclic dimensions.

# Output:

* `stepsize`: spatial (and temporal) shift in grid points between subdomains for every dimension (?)
* `overlapping`: number of overlapping grid points for every dimension
* `isdirect`: true is the direct solver is activated

"""
function DIVAnd_fittocpu(Lpmnrange,gridsize,latercsteps,moddim,MEMTOFIT)
    #################################################################################
    # Number of dimensions

    n = size(Lpmnrange,1)
    fudgefac=MEMTOFIT/16.

    if moddim==[]
        moddim=zeros(n)
    end

    # Some tweaking parameters #####

    # Fraction of domain size at which windowing is attempted if len is smaller
    lfactor=0.2

    # How wide is the overlap in terms of number of length scales
    factoroverlap=3.3

    biggestproblemitern=[500*500 500*500 50*50*50 170*170*6*12]*fudgefac
    biggestproblemdirectn=[200*200 200*200 50*50*40 50*50*10]*fudgefac

    biggestproblemiter=biggestproblemitern[min(n,4)]
    biggestproblemdirect=biggestproblemdirectn[min(n,4)]

    @debug "biggestproblemdirectn $biggestproblemdirectn"

    # Default
    stepsize=collect(gridsize)
    overlapping=zeros(Int,n)

    # Test if cutting is necessary:
    problemcut = prod(2*overlapping+stepsize)
    isdirect = problemcut < biggestproblemdirect

    @debug "problemcut: $problemcut"
    @debug "biggestproblemdirect: $biggestproblemdirect"

    if isdirect
        return stepsize,overlapping,isdirect
    end


    #####################################################################################
    # Define overlapping and stepsize

    # Now hardwired windowing over x,y
    # For the moment on z: overlapping of two layers and stepping of 4 keeps the size small enough
    # and overhead a factor of two

    # For time: no windowing for the moment neither

    biggestproblem= biggestproblemiter

    higherdims=1

    if n>2
        stepsize[3]=2;
        overlapping[3]=2;
        higherdims=prod(stepsize[3:end]+2*overlapping[3:end])
    end

    biggestproblem=biggestproblem/higherdims


    # problemsize is the number additional grid point appended to
    # a subdomain to make the domains overlap
    problemsize=1

    # nwd calculates the number of dimensions on which tiling/windowing is done
    nwd=0
    for i=1:min(n,2)
        # if length scale is small compared to domain size
        if Lpmnrange[i][2]<   lfactor*gridsize[i]
            if moddim[i]==0
                overlapping[i]=Int(ceil( factoroverlap*Lpmnrange[i][2]   ))
                problemsize=problemsize*overlapping[i]


                nwd=nwd+1
            else
                problemsize=problemsize*gridsize[i]
            end
        else
            problemsize=problemsize*gridsize[i]
        end

        #
    end

    #problemsize=problemsize/prod(latercsteps[1:2])
	
    # Take into account overhead due to multiple storage
    problemsize=problemsize/sqrt(prod(latercsteps[1:2]))    
    epsilon = 1E-6

    # tries to get the maximum multiplication factor with respect to the overlapping which can be applied 
    # to get the actual useful window size excluding the overlapping.
    if nwd>0
        epsilon=(float(biggestproblem)/float(problemsize))^(1.0/nwd)-2.0
    end

    if epsilon<=0
	if nwd>0
        @warn "Problem size probably too big for the memory defined " *
            "(epsilon_fittocpu = $epsilon, problemsize = $problemsize, nwd = $nwd" *
            "overlapping = $overlapping). Will try to continue anyway. " *
            "Consider to increase MEMTOFIT."
	end
        epsilon=1E-6
    end

    for i=1:min(n,2)
        # if length scale is small compared to domain size
        if Lpmnrange[i][2]<   lfactor*gridsize[i]
            if moddim[i]==0
                stepsize[i]=Int(ceil( epsilon*factoroverlap*Lpmnrange[i][2]))
                # Limit stepsize, if limited force zero overlap
                stepsize[i] = min(stepsize[i],gridsize[i])
		if gridsize[i]==stepsize[i]
			overlapping[i]=0
		end
			
            end
        end
    end

    # Limit window size to gridsize
    winsize = min.(2*overlapping+stepsize,gridsize)
    @debug "winsize: $winsize"

    doesitfit = prod(winsize) < biggestproblemiter

    # Before returning, check if by chance the windows are now small enough to even allow for a direct solver
    isdirect = prod(winsize) < biggestproblemdirect

    ####################################
    #Force direct solver if you want by uncommenting next line
    # isdirect = true

    return stepsize,overlapping,isdirect
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

# LocalWords:  fi DIVAnd pmn len diag CovarParam vel ceil moddim fracdim
