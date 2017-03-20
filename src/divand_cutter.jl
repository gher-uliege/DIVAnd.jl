"""


windowlist,csteps,lmask,alphapc = divand_cutter(Lpmnrange,gridsize,moddim=[]);

# Creates a list of windows for subsequent domain decomposition
# Also calculates already the subsampling steps csteps for the preconditionners
# as well as the mask lmask to apply to the length scales in the preconditionner, allowing to reduce
# the problem size

# Input:

* `Lpmnrange`:

* `gridsize`: number of points in each direction (size(mask))

* `moddim`:



# Output:

* `windowlist`: Array of tuples (iw1 iw2 ...)
* `csteps` : Array of steps for the coarse grid preconditionner
* `lmask` : Array of multiplication factors for length scale of preconditionner
* `alphapc` : Norm defining coefficients for preconditionner

"""


function divand_cutter(Lpmnrange,gridsize,moddim=[])




    # Some tweaking parameters #####

    # Minimum number of points per length scale for the preconditionner can be fractional
    minimumpointsperlpc=3.0






    #################################################################################
    # Number of dimensions
    n = size(Lpmnrange)[1]

    ######################################################
    # Calculate the sampling steps for the preconditionner
    #
    # If you want a direct solver put csteps zero
    ######################################################
    csteps=ones(Int,n);
    for i=1:n
        nsamp=Int(floor(Lpmnrange[i][1]/minimumpointsperlpc));
        if nsamp>1
            csteps[i]=minimum([nsamp,3])
        end
    end
    @show csteps
    #################################################################
    # Decide which directions are not coupled during preconditionning
    lmask=ones(n)

    # For the moment hardwired decoupling on z only


    lmask=ones(n)


    alphapc=[1 2 1]

    if n==4
        lmask[3]=0
        lmask=lmask*1.25715/0.69315
    end


    if n==3
        lmask=lmask*1.25715/0.69315
    end




    #####################################################################################
    # Define overlapping and stepsize


    # For the moment on z: overlapping of two layers and stepping of 4 keeps the size small enough
    # and overhead a factor of two

    # For time: if periodic, do windowing, otherwise as for x and y ?


    stepsize,overlapping,isdirect=divand_fittocpu(Lpmnrange,gridsize,csteps,moddim)

    if isdirect
        # Indiciate to the calling one that direct method can be used on windows
        csteps=0*csteps
        warn("Testing forced jog")
    end


    #####################################################################################
    # Now prepare the window infos

    # Keep pointers to place where solution is found in output grid

    # window corners in main problem
    iw1=zeros(Int,n)
    iw2=zeros(Int,n)

    # Solution indexes in submodel
    isol1=zeros(Int,n)
    isol2=zeros(Int,n)

    # Range for final storage
    istore1=zeros(Int,n)
    istore2=zeros(Int,n)

    # Running pointer
    ij=ones(Int,n)


    #
    # need to subtract two overlap regions from total size to determine the number of tiles
    #subsz = ([ceil(Int,sz[i] / stepsize[i]) for i = 1:ndims(mask)]...)

    subsz = ([ceil(Int,(gridsize[i]-2*overlapping[i]) / stepsize[i]) for i = 1:n]...)

    ntiles=prod(subsz)


    windowlist=Array(Tuple,ntiles)

    iw=0
    for cr in CartesianRange(subsz)
        ij = [[(cr[i]-1)*stepsize[i]+1  for i = 1:n]...]

        iw=iw+1
        #
        for nd=1:n

            iw1[nd]=ij[nd]

            # For normal tiles take middle part
            isol1[nd]=overlapping[nd]+1

            #JM
            # For first tile take up to the left limit
            if ij[nd]==1
                isol1[nd]=1
            end

            # Tile size is stepsize[nd]+2*overlapping[nd]
            iw2[nd]=ij[nd]+stepsize[nd]+2*overlapping[nd]-1;


            # For normal tiles take middel part
            isol2[nd]=overlapping[nd]+stepsize[nd]

            if iw2[nd]>=gridsize[nd]
                # Last box
                iw2[nd]=gridsize[nd]
                # Take up to the right border
                isol2[nd]=iw2[nd]-iw1[nd]+1
            end

            istore1[nd]=ij[nd];

            #JM
            if ij[nd]>1
                istore1[nd]=istore1[nd]+overlapping[nd]
            end

            istore2[nd]=istore1[nd]+isol2[nd]-isol1[nd];


        end
        # end loop over all dimensions

        # Need a deepcopy here otherwise last values everywhere as only memory adress would be used

        windowlist[iw]=deepcopy((iw1,iw2,isol1,isol2,istore1,istore2,))

    end
    # end loop over all windows








    return windowlist,csteps,lmask,alphapc



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
