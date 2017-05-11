"""
Compute a variational analysis of arbitrarily located observations.

fi,s = divandjog(mask,pmn,xi,x,f,len,epsilon2,csteps,lmask; alphapc=[1 2 1], otherargs...);

Perform an n-dimensional variational analysis of the observations `f` located at
the coordinates `x`. The array `fi` represent the interpolated field at the grid
defined by the coordinates `xi` and the scales factors `pmn`.

# Input:
* Same parameters as for divarun.
        * Two additional parameters:
                * csteps: array of ndims values defining the sampling steps for the preconditionner
                * lmask: array of ndims mutilplications factors for length scales
        * One additional optiional parameter
                * alphapc: The coefficients for the norm used in the preconditionner



# Output:
*  `fi`: the analysed field
*  `s`: structure with an array `s.P` representing the analysed error covariance

# Note:


"""



function divandjog(mask,pmn,xi,x,f,Labs,epsilon2,csteps,lmask; alphapc=[],otherargs...
                   )
    #

    n=ndims(mask)
    nsteps=csteps




    if sum(nsteps)==0
        #####################################################
        # Run normal direct version
        #####################################################
        finter,sinter=divandrun(mask,pmn,xi,x,f,Labs,epsilon2; otherargs...)
        return finter,sinter
    end


    ##########################################################################################
    # Capture here the case sum(csteps)==n
    if sum(nsteps)==n

        # In this case, can I just define HI=eye(), or something as sparse
        @show sum(nsteps)
        # HI=speye(prod(size(mask)))
        # But in reality not needed if tests on sum(ntests) done so maybe better really separate the two cases
        # Would also simplify the preconditionner calculations !

        # Normally the case with no subsampling is only usefull if alphapc is defined explicitely
        # or lmask decouples some directions
        # (and different values of alpha from the normally used one are defined but this is not tested here)
        @show(prod(lmask))
        if ((alphapc==[]) & (prod(lmask)>0))
            warn("divajog called with no coarsening and normal alpha")
            warn("pcg without preconditionner will be used")
            warn("to force use of direct solver put csteps to zero")
            fi,si=divandrun(mask,pmn,xi,x,f,Labs,epsilon2; otherargs...,inversion=:pcg)
            return fi,si
        else
            # Case where the grid is not subsampled but different norms used for preconditionning

            @show sum(nsteps)

            if isa(Labs,Tuple)
                Labsc=Labs
            else
                Labsc=(Labs*ones(n)...);
            end
            Labsccut=([Labsc[i]*lmask[i] for i=1:n]...)
            # Run model with simplified norm
            fc,sc=divandrun(mask,pmn,xi,x,f,Labsccut,epsilon2; otherargs...,alpha=alphapc)

            if !(sc==0)
                # TEST makes sure there are values in the coarse resolution solution
                # Preconditionner core
                scP=sc.P;
                figuess=fc
            else
                scP=1
                figuess=zeros(size(mask))
            end

            # Try to clean up some memory here
            s=0
            fc=0
            gc()

            # To compensate for the missing correlations in scP
            # tolerance on the gradient A x - b
            tol = 2e-3
            maxiter=10*Int(ceil(sqrt(prod(size(mask)))))
            pcargs = [(:tol, tol),(:maxit,maxiter)]

            diagshift=0.0001;

            # Preconditionner function
            function compPCa(iB,H,R)
                function fun!(x,fx)
                    fx[:] = diagshift*x+scP*x;
                end
                return fun!
            end

            # Then run with normal resolution and preconditionner

            fi,si=divandrun(mask,pmn,xi,x,f,Labs,epsilon2; otherargs...,pcargs...,inversion=:pcg,operatortype=Val{:MatFun},compPC = compPCa, fi0 =figuess)
            #fi,si=divandrun(mask,pmn,xi,x,f,Labs,epsilon2; otherargs...,pcargs...,inversion=:pcg,compPC = compPCa, fi0 =figuess)
            return fi,si
        end
    end


    #########################################################################################
    # Now the more complicated case

    # Need to check for cyclic boundaries for later calculation of HI

    moddim=zeros(n);
    kwargs_dict = Dict(otherargs)
    if haskey(kwargs_dict,:moddim)
        moddim=kwargs_dict[:moddim]
    end
    iscyclic = moddim .> 0


    if sum(nsteps)>n


        ####### use preconditionner method with coarsening.




        #######################################
        # HOW TO TREAT COARSE GRID NOT COVERING FINE GRID ?
        # Artificially add of coordinates of last points in fine grid if the last step is beyond.
        # Slight inconsistency here for the error field as pmn was not adapted for last two points
        #######################################

        #coarsegridpoints=([1:nsteps[i]:size(mask)[i] for i in 1:n]...);
        #([unique(push!(collect(1:3:13),13)) for i in 1:4]...)

        # coarsegridpoints=([unique(push!(collect(1:nsteps[i]:size(mask)[i]),size(mask)[i])) for i in 1:n]...);

        # Test for better convergence making sure the expanded points are take in the preconditionner but then expande the grid only with alpha=0.5 since factor 2 typically applied later ?
        # Alphabc should be an array to do this properly ...



        coarsegridpoints=([sort(unique(push!(collect(2:nsteps[i]:size(mask)[i]-1),size(mask)[i],size(mask)[i]-1,1))) for i in 1:n]...);

        # If last point not reached add last point and just forget about incorrect metric  there ?


        #

        xic=([ x[coarsegridpoints...] for x in xi ]...);


        # Create a slightly expanded sea mask if nsteps are not 1


        maskf=deepcopy(mask)
        for ii=1:n
            if nsteps[ii]>1
                maskf=mapslices(dvmaskexpand,maskf,ii)
            end
        end
        maskc=maskf[coarsegridpoints...];



        # Now scale pmn by the step factors

        pmnc=([ (1.0/nsteps[i])*pmn[i][coarsegridpoints...] for i=1:length(pmn) ]...)

        # Check if Labs is a tuple of tuple; in this case also subsample

        if isa(Labs,Tuple)
            if isa(Labs[1],Number)
                Labsc=Labs;
            else
                Labsc=([ x[coarsegridpoints...] for x in Labs ]...);
            end
        else
            # Create a tuple of L for the coarse grid; needed to be able to put some of them to zero
            #  Labsc=Labs;
            Labsc=(Labs*ones(n)...);
        end


        # Create HI only of subsampled, otherwise use eye
        # Now prepare HI do go from the coarse grid to the fine grid. To do so
        # interprete de fine grid coordinates as those of pseudo-obs and use divandtoos
        # Need the statevector strucure for the fine grid to go from grid to array

        svf = statevector_init((mask,))

        # For each coordinate in the tuplet xi, go from grid representation to tuplet
        # to have the pseudo-data coordinates


        Ic = localize_separable_grid(([statevector_pack(svf,(x,)) for x in xi]...),maskc,xic);

        # Create fractional indexes of these data points in the coarse grid

        HI,outc,outbboxc = sparse_interp(maskc,Ic,iscyclic);
        HI = HI * sparse_pack(maskc)';

        ####################################
        # Advection constraint extracted
        ####################################



        # Search for velocity argument:
        jfound=0
        for j=1:size(otherargs)[1]
            if otherargs[j][1]==:velocity
                jfound=j
                break
            end
        end

        if jfound>0
            # modify the parameter only in the coarse model
            otherargsc=deepcopy(otherargs)
            otherargsc[jfound]=(:velocity,([ x[coarsegridpoints...] for x in otherargs[jfound][2] ]...))
        else
            otherargsc=otherargs
        end



        # For other constraints:
        # TODO TODO TODO TODO IF divandjog should accept other constraints
        # Here should be straightfoward replace C by C*HI on the constraint structure





        # Prepare run of the coarse grid problem





        # For the coarse model, slightly adapth alphabc assuming a typical ratio of 4 is used
        # Search for alphabc argument:
        kfound=0
        for j=1:size(otherargs)[1]
            if otherargs[j][1]==:alphabc
                kfound=j
                break
            end
        end
        if kfound>0
            if jfound==0
                otherargsc=deepcopy(otherargs)
            end
            # modify the parameter only in the coarse model
            otherargsc[kfound]=(:alphabc,0.25)
        else
            #       warn("Need to expand")
            otherargsc=vcat(otherargsc,(:alphabc,0.25))
        end
        @show otherargsc
        @show alphapc

        # maybe try another norm using btrunc=3 or even btrunc=2 but full L here for the coarser version ?

        # Preconditionner with desactivated correlations in some directions

        # Hardcoded new test


        Labsccut=([Labsc[i]*lmask[i] for i=1:n]...)

        Labsccut=Labsc


        # Run coarse resolution model
        #               fc,sc=divandrun(maskc,pmnc,xic,x,f,Labsccut,epsilon2; otherargsc...,alpha=alphapc,btrunc=2)
        fc,sc=divandrun(maskc,pmnc,xic,x,f,Labsccut,epsilon2; otherargsc...,btrunc=3)





        if !(sc==0)
            # TEST makes sure there are values in the coarse resolution solution
            # Take the fc coarse solution, pack to to statevector form
            # using sc.sv
            # Apply HI; this vector can also be used as a first guess for the PC

            # Preconditionner core
            scP=sc.P;

            xguess=HI*statevector_pack(sc.sv,(fc,));
            figuess,=statevector_unpack(svf,xguess)
        else

            # Only sea points in coarse grid
            scP=1
            figuess=zeros(size(mask))
        end

# Try to clean up some memory here
s=0
fc=0
gc()



# tolerance on the gradient A x - b
tol = 2e-3

maxiter=10*Int(ceil(sqrt(size(HI)[1])))
maxiter=1000
pcargs = [(:tol, tol),(:maxit,maxiter)]

# To compensate for the missing correlations in HI*scP*HI'

diagshift=0.006*(sqrt(size(HI)[1]/size(HI)[2])-1);
diagshift=0.03*(sqrt(size(HI)[1]/size(HI)[2])-1);
diagshift=0.04*(sqrt(size(HI)[1]/size(HI)[2])-1);
#diagshift=0.01
#                Z=randn(size(HI)[2],5);
#                diagshift=mean(diagMtCM(scP,Z)./diag(Z'*Z))
#        diagshift=0.021*diagshift

# Preconditionner function
function compPC(iB,H,R)
    function fun!(x,fx)
        fx[:] = diagshift*x+HI*(scP*(HI'*x))
    end
    return fun!
end





# Then run with normal resolution and preconditionner

#fi,si=divandrun(mask,pmn,xi,x,f,Labs,epsilon2; otherargs...,pcargs...,inversion=:pcg,operatortype=Val{:MatFun},compPC = compPC, fi0 =figuess)
fi,si=divandrun(mask,pmn,xi,x,f,Labs,epsilon2; otherargs...,pcargs...,inversion=:pcg,compPC = compPC, fi0 =figuess,btrunc=2)


# Some ideas for error calculations
#fi,si=divandrun(mask,pmn,xi,x,f,Labs,epsilon2; otherargs...,pcargs...,inversion=:pcg,compPC = divand_pc_sqrtiB, fi0 =xguess)
#errfield=diagMtCM(sc.P,HI')
#erri,=statevector_unpack(si.sv,errfield)
# For error field based on coarse one, use divand_filter3 with ntimes=Int(ceil(mean(nsteps)))
# First test

# Possible optimization for a climatology production, for each tile return diagshift and niter with some random changes in diagshift and search
# for optimum. In particular if several climatology runs are produced (slight changed parameters), the preliminary runs, without error calculations for example
# could provide good estimates ?


return fi,si
end
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
