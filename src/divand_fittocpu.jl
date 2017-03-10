"""


stepsize,overlapping,isdirect = divand_fittocpu(Lpmnrange,gridsize,moddim=[]);

# Creates a list of windows for subsequent domain decomposition
# Also calculates already the subsampling steps csteps for the preconditionners

# Input:

* `Lpmnrange`:

* `gridsize`: number of points in each direction (size(mask))

* `moddim`:



# Output:


"""


function divand_fittocpu(Lpmnrange,gridsize,moddim=[])


    #################################################################################
    # Number of dimensions
    n = size(Lpmnrange)[1]

    if moddim==[]
        moddim=zeros[n]
    end

    # Some tweaking parameters #####


    # Fraction of domain size at which windowing is attempted if len is smaller
    lfactor=0.2

    # How wide is the overlap in terms of number of length scales
    factoroverlap=3.3

	biggestproblemitern=[500*500 500*500 50*50*50 55*55*10*12]
	biggestproblemdirectn=[200*200 200*200 50*50*20 50*50*10]
	
	biggestproblemiter=biggestproblemitern[minimum([n,4])]
	biggestproblemdirect=biggestproblemdirectn[minimum([n,4])]
	
    #if n<3
    #    biggestproblemiter=500*500
    #    biggestproblemdirect=200*200
    #end
    #if n==3
    #    biggestproblemiter=50*50*50
    #    biggestproblemdirect=50*50*20
    #end
    #if n>3
    #    biggestproblemiter=55*55*10*12
    #    biggestproblemdirect=50*50*10
    #end








    # Default
    stepsize=collect(gridsize)
    overlapping=zeros(Int,n)

    # Test if cutting is necessary:
    isdirect=(prod(2*overlapping+stepsize)<biggestproblemdirect)

    if isdirect
        return stepsize,overlapping,isdirect
    end




    # Unfortunataly for the moment the problem is memory bound by the unsampled grid.

    laterscales=ones(n)

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

    problemsize=1


    nwd=0
    for i=1:minimum([n,2])


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

    problemsize=problemsize/prod(laterscales[1:2])


    if nwd>0
        epsilon=(float(biggestproblem)/float(problemsize))^(1.0/nwd)-2.0
    end
    if epsilon<0
        warn("SO what $epsilon $problemsize $nwd $overlapping")
        epsilon=1E-6
    end

    for i=1:minimum([n,2])
        # if length scale is small compared to domain size
        if Lpmnrange[i][2]<   lfactor*gridsize[i]
            if moddim[i]==0
                stepsize[i]=Int(ceil( epsilon*factoroverlap*Lpmnrange[i][2]   ))
            end
        end
    end




    doesitfit=(prod(2*overlapping+stepsize)<biggestproblemiter)





    # Before returning, check if by chance the windows are now small enough to even allow for a direct solver

    isdirect=(prod(2*overlapping+stepsize)<biggestproblemdirect)

    ####################################
    #Force direct solver if you want by uncommenting next line
    # isdirect=(0<1)



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

# LocalWords:  fi divand pmn len diag CovarParam vel ceil moddim fracdim
