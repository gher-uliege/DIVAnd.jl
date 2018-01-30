"""
Compute a variational analysis of arbitrarily located observations.

fi,s = divandgo(mask,pmn,xi,x,f,len,epsilon2,errormethod; ...);

Perform an n-dimensional variational analysis of the observations `f` located at
the coordinates `x`. The array `fi` represent the interpolated field at the grid
defined by the coordinates `xi` and the scales factors `pmn`.

# Input:

As for divandrun but as a higher level routine which will automatically create windowing etc
it also include the definition of the errormethod

* errormethod : :cpme (clever poormans method), :none or :exact


# Output:
*  `fi`: the analysed field
*  `s`: structure with an array `s.P` representing the analysed error covariance

"""



function divandgo(mask,pmn,xi,x,f,Labs,epsilon2,errormethod=:cpme; otherargs...
                  )







    n=ndims(mask)
    # Needed to make sure results are saved.
    fi=0
    s=0

    # Need to check for cyclic boundaries

    moddim=zeros(n);

    kwargs_dict = Dict(otherargs)

    if haskey(kwargs_dict,:moddim)
        moddim=kwargs_dict[:moddim]
    end
	# optional argument for window fitting.
	MEMTOFIT=16
	if haskey(kwargs_dict,:MEMTOFIT)
        MEMTOFIT=kwargs_dict[:MEMTOFIT]
    end

    # DOES NOT YET WORK WITH PERIODIC DOMAINS OTHER THAN TO MAKE SURE THE DOMAIN IS NOT CUT
    # IN THIS DIRECTION. If adapation is done make sure the new moddim is passed to divandrun
    # General approach in this future case prepare window indexes just taking any range including negative value and
    # apply a mod(myindexes-1,size(mask)[i])+1 in direction i when extracting
    # for coordinates tuples of the grid (xin,yin, .. )  and data (x,y)
    # in the direction, shift coordinates and apply modulo mod(x-x0+L/2,L)
    #




    # Analyse rations l/dx etc

    Lpmnrange = divand_Lpmnrange(pmn,Labs)

    # Create list of windows, steps for the coarsening during preconditioning and mask for lengthscales to decoupled directions during preconditioning
	@show moddim,MEMTOFIT
	
    windowlist,csteps,lmask,alphanormpc = divand_cutter(Lpmnrange,size(mask),moddim,MEMTOFIT)

	@show size(mask),size(windowlist)

    # For parallel version declare SharedArray(Float,size(mask)) instead of zeros() ? ? and add a @sync @parallel in front of the for loop ?
    # Seems to work with an addprocs(2); @everywhere using divand to start the main program. To save space use Float32 ?
    #fi=zeros(size(mask));

    #fi=SharedArray(Float64,size(mask));
    #erri=SharedArray(Float64,size(mask));
    fi=SharedArray{Float32}(size(mask));
	if errormethod==:none
        erri=eye(1)
        else
	    erri=SharedArray{Float32}(size(mask));
    end

    @sync @parallel for iwin=1:size(windowlist)[1]

        iw1=windowlist[iwin][1]
        iw2=windowlist[iwin][2]
        isol1=windowlist[iwin][3]
        isol2=windowlist[iwin][4]
        istore1=windowlist[iwin][5]
        istore2=windowlist[iwin][6]


		windowpointssol=([isol1[i]:isol2[i] for i in 1:n]...);
		windowpointsstore=([istore1[i]:istore2[i] for i in 1:n]...);



        warn("Test window $iw1 $iw2 $isol1 $isol2 $istore1 $istore2 ")


        windowpoints=([iw1[i]:iw2[i] for i in 1:n]...);

        @show size(windowpoints[1])

        #################################################
        # Need to check how to work with aditional constraints...
        #################################################

        #################################
        # Search for velocity argument:
        jfound=0
        for j=1:size(otherargs)[1]
            if otherargs[j][1]==:velocity
                jfound=j
                break
            end
        end






        if jfound>0
            warn("There is an advection constraint; make sure the window sizes are large enough for the increased correlation length")
            # modify the parameter
            otherargsw=deepcopy(otherargs)
            otherargsw[jfound]=(:velocity,([ x[windowpoints...] for x in otherargs[jfound][2] ]...))
        else
            otherargsw=otherargs
        end







        # If C is square then maybe just take the sub-square corresponding to the part taken from x hoping the constraint is a local one ?
        #





        # If C projects x on a low dimensional vector: maybe C'C x-C'd as a constraint, then pseudo inverse and woodbury to transform into a similar constraint but on each subdomain
        # Would for example replace a global average constraint to be replaced by the same constraint applied to each subdomain. Not exact but not too bad neither






        fw=0
        s=0

        xiw=deepcopy(([ x[windowpoints...] for x in xi ]...))


        # NEED TO CATCH IF Labs is a tuple of grid values; if so need to extract part of interest...

        Labsw=Labs
        if !isa(Labs,Number)
            if !isa(Labs[1],Number)
                Labsw= ([ x[windowpoints...] for x in Labs ]...)
            end
        end

        kfound=0
        for j=1:size(otherargs)[1]
            if otherargs[j][1]==:alphabc
                kfound=j
                break
            end
        end

        if 3==2
            if kfound>0
                # modify the parameter only in the window model
                if jfound>0
                    # Ok this is alreay a copy and you can change it
                    otherargsw[kfound]=(:alphabc,1)
                else
                    otherargsw=deepcopy(otherargs)
                    otherargsw[kfound]=(:alphabc,1)
                end
            else
                warn("Need to expand")
                otherargsw=vcat(otherargsw,(:alphabc,1))
            end
        end

        # If you want to change another alphabc, make sure to replace it in the arguments, not adding them since it already might have a value
        # Verify if a direct solver was requested from the demain decomposer
        if sum(csteps)>0
		    fw,s=divandjog(mask[windowpoints...],([ x[windowpoints...] for x in pmn ]...),xiw,x,f,Labsw,epsilon2,csteps,lmask;alphapc=alphanormpc, otherargsw... )
            fi[windowpointsstore...]= fw[windowpointssol...];
			if errormethod==:cpme
			    fw=0
                s=0
                gc()
                # Possible optimization here: use normal cpme (without steps argument but with preconditionner from previous case)
                errw=divand_cpme(mask[windowpoints...],([ x[windowpoints...] for x in pmn ]...),xiw,x,f,Labsw,epsilon2;csteps=csteps,lmask=lmask,alphapc=alphanormpc, otherargsw... )
            end
            # for errors here maybe add a parameter to divandjog ? at least for "exact error" should be possible; and cpme directly reprogrammed here as well as aexerr ? assuming s.P can be calculated ?
        else
            # Here would be a natural place to test which error fields are demanded and add calls if the direct method is selected
            fw,s=divandrun(mask[windowpoints...],([ x[windowpoints...] for x in pmn ]...),xiw,x,f,Labsw,epsilon2; otherargsw...)
			fi[windowpointsstore...]= fw[windowpointssol...];
            if errormethod==:cpme
                errw=divand_cpme(mask[windowpoints...],([ x[windowpoints...] for x in pmn ]...),xiw,x,f,Labsw,epsilon2; otherargsw... )
            end
        end

        

        # Cpme: just run and take out same window

        # AEXERR: just run and() take out same window

        if errormethod==:exact
            # EXERR: P only to points on the inner grid, not the overlapping one !
            # Initialize errw everywhere,
            errw=0.*fw
            # For packing, take the stavevector returned except when sum(csteps)>n
            # in this case recreate a statevector
            if sum(csteps)==n
                sverr=statevector_init((mask[windowpoints...],))
            else
                sverr=s.sv
            end
            svn=size(sverr)[1]
            @show size(sverr)
            errv=statevector_pack(sverr,(errw,))
            # Loop over window points. From grid index to statevector index so that ve is
            # zero exect one at that index. Then calculate the error and store it in the the
            # sv representation of the error

            @show windowpoints[1]
            for gridindex in windowpoints
                @show gridindex
                ei=zeros(svn)
                ind = statevector_sub2ind(svn,gridindex)
                ei[ind]=1
                #  HIP=HI'*ei
                #  errv[ind]=diagMtCM(sc.P,HIP)
            end
            errw=statevector_unpack(svn,errv)
            # Maybe better fill in first HIP = HI'*[ ... ei ...]
            # then something as errfield=diagMtCM(sc.P,HIP)
            # at the end of the loop, unpack the sv error field into errw
            # End error fields
        end


# copy, deepcopy or just = ???




		if errormethod==:none
            erri=eye(1);
			else
			erri[windowpointsstore...]=errw[windowpointssol...];
		end

end

# When finished apply an nd filtering to smooth possible edges, particularly in error fields.
# it also makes the shared array possible to save in netCDF??
fi=divand_filter3(fi,NaN,2)
erri=divand_filter3(erri,NaN,3)


return fi,erri



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
