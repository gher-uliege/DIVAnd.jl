"""
Compute a variational analysis of arbitrarily located observations.

fi,s = DIVAndjog(mask,pmn,xi,x,f,len,epsilon2,csteps,lmask; alphapc=[1,2,1], otherargs...);

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


"""
function DIVAndjog(mask,pmn,xi,x,f,Labs,epsilon2,csteps,lmask,pcmethod=1; alphapc=[],otherargs...
                   )
    #

    n=ndims(mask)
    nsteps=csteps




    if sum(nsteps)==0
        #####################################################
        # Run normal direct version
        #####################################################
        finter,sinter=DIVAndrun(mask,pmn,xi,x,f,Labs,epsilon2; otherargs...)
        return finter,sinter
    end

    ##########################################################################################
    # Capture here the case sum(csteps)==n
    if sum(nsteps)==n

        if ( (pcmethod==0) & ((alphapc==[]) & (prod(lmask)>0)) )
            @warn "divajog called with no coarsening and no alphapc defined"
            #@warn "pcg without preconditionner will be used"
            #@warn "to force use of direct solver put csteps to zero"
		    tol = 2e-3
			maxiter=100*Int(ceil(sqrt(prod(size(mask)))))

            fi,si=DIVAndrun(mask,pmn,xi,x,f,Labs,epsilon2; otherargs...,
                            tol = tol,
                            maxit = maxiter,
                            inversion=:pcg)
            return fi,si

        else


			tol = 2e-3
			maxiter=10*Int(ceil(sqrt(prod(size(mask)))))
			diagshift=0.0000001;

            # Case where the grid is not subsampled but different norms used for preconditionning



            if isa(Labs,Tuple)
                Labsc=Labs
            else
                Labsc=(Labs*ones(n)...,);
            end

			methodpc=pcmethod
            fi=0
			si=0

			if methodpc==1
			    diagshift=0.0001;
				Labsccut=([Labsc[i]*lmask[i] for i=1:n]...,)
				# Run model with simplified norm
				fc,sc=DIVAndrun(mask,pmn,xi,x,f,Labsccut,epsilon2; otherargs...,alpha=alphapc)
                scP=1
				figuess=zeros(size(mask))
				xguess=1
				if !(sc==0)
                    # TEST makes sure there are values in the coarse resolution solution
                    # Preconditionner core
					scP=sc.P;
					figuess=fc
				end

				# Try to clean up some memory here
				sc=0
				fc=0
				#gc()

				# To compensate for the missing correlations in scP
				# tolerance on the gradient A x - b


				# Preconditionner function
                function compPCa(iB,H,R)
                    function fun!(x,fx)
						#fx[:] = diagshift*x+scP*x;
						fx[:]=scP*x
						fx[:]=BLAS.axpy!(diagshift,x,fx)

                    end
                    return fun!
                end
				# Then run with normal resolution and preconditionner

				fi,si=DIVAndrun(mask,pmn,xi,x,f,Labs,epsilon2; otherargs...,
                                tol = tol,
                                maxit = maxiter,
                                inversion=:pcg,
                                compPC = compPCa, fi0 =figuess,btrunc=2)

			end


			if methodpc == 2
				lmask1 = 0 .* lmask;
				lmask1[1:2] .= 1.0;
				Labsccut=([Labsc[i]*lmask1[i] for i=1:n]...,)
				fc,sc=DIVAndrun(mask,pmn,xi,x,f,Labsccut,epsilon2; otherargs...)
				PC2=1
				xguess=1
				if !(sc==0)
					# TEST makes sure there are values in the coarse resolution solution
					# Take the fc coarse solution, pack to to statevector form
					# using sc.sv
					# Apply HI; this vector can also be used as a first guess for the PC

					# Preconditionner core
					PC2=sc.P;
					xguess=statevector_pack(sc.sv,(fc,));
				end
				# Try to clean up some memory here
				sc=0
				fc=0
				#gc()
				lmask1=0.0 .* lmask;
				lmask1[3:end] .= 1.0./1.42;

				Labsccut=([Labsc[i]*lmask1[i] for i=1:n]...,)
				PC1=1
				sc=0
				if size(lmask)[2]>2
					fc,sc=DIVAndrun(mask,pmn,xi,x,f,Labsccut,10000.; otherargs...)
				end
				if !(sc==0)
				    # TEST makes sure there are values in the coarse resolution solution
				    # Take the fc coarse solution, pack to to statevector form
				    # using sc.sv
				    # Apply HI; this vector can also be used as a first guess for the PC

				    # Preconditionner core
					PC1=sc.P;
				end
			    # Try to clean up some memory here
				sc=0
				fc=0
				#gc()




			    # tolerance on the gradient A x - b





			    xguess=PC1*(PC1*xguess);
			    xr=randn(size(xguess)[1],1)
			    scalef=(xr'*(PC1*(PC2*(PC1*xr))))./(xr'*xr)
			    scalefter=(xr'*(PC2*xr))./(xr'*xr)



			    scalef2=scalefter[1]/scalef[1]


			    xguess=xguess*scalef2
			    svf = statevector_init((mask,))
			    figuess,=statevector_unpack(svf,xguess)

			    function compPC4(iB,H,R)
                    function fun!(x,fx)
                        fx[:] = diagshift*x+scalef2*(PC1*(PC2*(PC1*x)))
                    end
                    return fun!
			    end

			    fi,si=DIVAndrun(mask,pmn,xi,x,f,Labs,epsilon2; otherargs...,
                                tol = tol,
                                maxit = maxiter,
                                inversion=:pcg,compPC = compPC4, fi0 =figuess,btrunc=2)


			end



			if methodpc==3
			    # same idea is for 3 but instead of trying to find an L such that B2 is B use directly decomposition of B!
				lmask1=0.0.*lmask;
				lmask1[1:2] .= 1.0;
				Labsccut=([Labsc[i]*lmask1[i] for i=1:n]...,)
				fc,sc=DIVAndrun(mask,pmn,xi,x,f,Labsccut,epsilon2; otherargs...)
				PC2=1
				xguess=1
				if !(sc==0)
					# TEST makes sure there are values in the coarse resolution solution
					# Take the fc coarse solution, pack to to statevector form
					# using sc.sv
					# Apply HI; this vector can also be used as a first guess for the PC

					# Preconditionner core
					PC2=sc.P;
					xguess=statevector_pack(sc.sv,(fc,));
				end
				# Try to clean up some memory here
				sc=0
				fc=0
				#gc()
				lmask1=0.0.*lmask;
				lmask1[3:end] .= 1.0;
				Labsccut=([Labsc[i]*lmask1[i] for i=1:n]...,)
				PC1=1
				sc=0

				if size(lmask)[2]>2
					fc,sc=DIVAndrun(mask,pmn,xi,x,f,Labsccut,10000.; otherargs...)
				end

				if !(sc==0)
				    # TEST makes sure there are values in the coarse resolution solution
				    # Take the fc coarse solution, pack to to statevector form
				    # using sc.sv
				    # Apply HI; this vector can also be used as a first guess for the PC

				    # Preconditionner core
					PC1=sc.P;
				end
			    # Try to clean up some memory here
				sc=0
				fc=0
				#gc()

			    # tolerance on the gradient A x - b

			    xguess=(PC1*xguess);
			    xr=randn(size(xguess)[1],1)
			    newsc=[1]
			    #scalef=(xr'*(PC1.factors[:UP]\(PC2*(PC1.factors[:PtL]\xr))))./(xr'*xr)
			    if PC1==1
			        scalef=(xr'*((PC2*(xr))))./(xr'*xr)
			    else
                    #			scalef=(xr'*(PC1.factors[:PtL]\(PC2*(PC1.factors[:UP]\xr))))./(xr'*xr)
  			        scalef=(xr'*(PC1.factors[:UP]\(PC2*(PC1.factors[:PtL]\xr))))./(xr'*xr)
			        newsc=(xr'*(PC1.factors[:UP]\((PC1.factors[:PtL]\xr))))./(xr'*xr)

			    end
			    scalefter=(xr'*(PC2*xr))./(xr'*xr)



			    scalef2=scalefter[1]/scalef[1]

			    scalef2=1.0./newsc[1]

			    xguess=xguess*scalef2
			    svf = statevector_init((mask,))
			    figuess,=statevector_unpack(svf,xguess)

			    function compPC4bis(iB,H,R)
                    function fun!(x,fx)
                        #fx[:] = diagshift*x+scalef2*(PC1.factors[:PtL]\(PC2*(PC1.factors[:UP]\x)))
                        fx[:] = diagshift*x+scalef2*(PC1.factors[:UP]\(PC2*(PC1.factors[:PtL]\x)))
                    end
                    return fun!
			    end
			    function compPC4bisb(iB,H,R)
                    function fun!(x,fx)
                        fx[:] = diagshift*x+scalef2*((PC2*(x)))
                    end
                    return fun!
			    end
			    if PC1==1
			        fi,si=DIVAndrun(
                        mask,pmn,xi,x,f,Labs,epsilon2; otherargs...,
                        tol = tol,
                        maxit = maxiter,
                        inversion=:pcg,compPC = compPC4bisb, fi0 =figuess,btrunc=2)
			    else
			        fi,si=DIVAndrun(
                        mask,pmn,xi,x,f,Labs,epsilon2; otherargs...,
                        tol = tol,
                        maxit = maxiter,
                        inversion=:pcg,compPC = compPC4bis, fi0 =figuess,btrunc=2)
			    end

			end




			if methodpc==4
			    # same idea is for 3 but instead of trying to find an L such that B2 is B use directly decomposition of B!
				lmask1=0.0.*lmask;
				lmask1[1:2] .= 1.0;
				Labsccut=([Labsc[i]*lmask1[i] for i=1:n]...,)
				fc,sc=DIVAndrun(mask,pmn,xi,x,f,Labsccut,epsilon2; otherargs...)
				PC2=1
				xguess=1
				if !(sc==0)
					# TEST makes sure there are values in the coarse resolution solution
					# Take the fc coarse solution, pack to to statevector form
					# using sc.sv
					# Apply HI; this vector can also be used as a first guess for the PC

					# Preconditionner core
					PC2=sc.P;
					xguess=statevector_pack(sc.sv,(fc,));
				end
				# Try to clean up some memory here
				sc=0
				fc=0
				#gc()
				lmask1=0.0.*lmask;
				lmask1[3:end] .= 1.0;
				Labsccut=([Labsc[i]*lmask1[i] for i=1:n]...,)
				# Try to get iB by using iterative solved stopped at one
				maxiterb=1

				PC1=1
				sc=0

				if size(lmask)[2]>2
					fc,sc=DIVAndrun(
                        mask,pmn,xi,x,f,Labsccut,10000.; otherargs...,
                        tol = tol,
                        maxit = maxiterb,
                        inversion=:pcg)
				end
				PC1=1
				if !(sc==0)
				    # TEST makes sure there are values in the coarse resolution solution
				    # Take the fc coarse solution, pack to to statevector form
				    # using sc.sv
				    # Apply HI; this vector can also be used as a first guess for the PC

				    # Preconditionner core
					PC1=sc.iB;
				end
			    # Try to clean up some memory here
				sc=0
				fc=0
				#gc()




			    # tolerance on the gradient A x - b





			    xr=randn(size(xguess)[1],1)
			    scalef=0.5*((xr'*(PC1*(PC2*xr)))+(xr'*(PC2*(PC1*xr))))./(xr'*xr)
			    #scalefter=(xr'*(PC2*xr))./(xr'*xr)



			    scalef2=0.03/(1+1.001*scalef[1])


			    svf = statevector_init((mask,))
			    figuess,=statevector_unpack(svf,xguess)

			    function compPC4ter(iB,H,R)
                    function fun!(x,fx)
                        fx[:] = diagshift*x+PC2*x-scalef2*(PC2*(PC1*(PC2*x)))
                    end
                    return fun!
			    end

			    fi,si=DIVAndrun(mask,pmn,xi,x,f,Labs,epsilon2; otherargs...,
                                tol = tol,
                                maxit = maxiter,
                                inversion=:pcg,compPC = compPC4ter, fi0 =figuess,btrunc=2)


			end

			if methodpc==5
			    # same idea is for 3 but instead of trying to find an L such that B2 is B use directly decomposition of B!
				lmask1=0.0.*lmask;
				lmask1[1:2] .= 1.0;
				Labsccut=([Labsc[i]*lmask1[i] for i=1:n]...,)
				fc,sc=DIVAndrun(mask,pmn,xi,x,f,Labsccut,epsilon2; otherargs...)
				PC2=1
				xguess=1
				if !(sc==0)
					# TEST makes sure there are values in the coarse resolution solution
					# Take the fc coarse solution, pack to to statevector form
					# using sc.sv
					# Apply HI; this vector can also be used as a first guess for the PC

					# Preconditionner core
					PC2=sc.P;
					xguess=statevector_pack(sc.sv,(fc,));
				end
				# Try to clean up some memory here
				sc=0
				fc=0
				#gc()




			    # tolerance on the gradient A x - b





			    scalef2=1.013

			    xguess=scalef2*xguess
			    svf = statevector_init((mask,))
			    figuess,=statevector_unpack(svf,xguess)

			    function compPC4quad(iB,H,R)
                    function fun!(x,fx)
                        fx[:] = diagshift*x+scalef2*(PC2*x)
                    end
                    return fun!
			    end

			    fi,si=DIVAndrun(mask,pmn,xi,x,f,Labs,epsilon2; otherargs...,
                                tol = tol,
                                maxit = maxiter,
                                inversion=:pcg,compPC = compPC4quad, fi0 =figuess,btrunc=2)


			end


			#			fs=statevector_pack(si.sv,(fi,));
			#			fgs=statevector_pack(si.sv,(figuess,));
			#			al=dot(fs,fgs)/dot(fgs,fgs)
			#			@show al,dot(fgs-fs,fgs-fs)/dot(fs,fs),dot(al*fgs-fs,al*fgs-fs)/dot(fs,fs)

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

        #coarsegridpoints=([1:nsteps[i]:size(mask)[i] for i in 1:n]...,);
        #([unique(push!(collect(1:3:13),13)) for i in 1:4]...,)

        # coarsegridpoints=([unique(push!(collect(1:nsteps[i]:size(mask)[i]),size(mask)[i])) for i in 1:n]...,);

        # Test for better convergence making sure the expanded points are take in the preconditionner but then expande the grid only with alpha=0.5 since factor 2 typically applied later ?
        # Alphabc should be an array to do this properly ...



        coarsegridpoints=([sort(unique(push!(collect(2:nsteps[i]:size(mask)[i]-1),size(mask)[i],size(mask)[i]-1,1))) for i in 1:n]...,);

        # If last point not reached add last point and just forget about incorrect metric  there ?


        #

        xic=([ x[coarsegridpoints...] for x in xi ]...,);


        # Create a slightly expanded sea mask if nsteps are not 1


        maskf=deepcopy(mask)
        for ii=1:n
            if nsteps[ii]>1
                maskf=mapslices(dvmaskexpand,maskf,dims = ii)
            end
        end
		# One further expansion, but maybe check if updates on the fly are really the way to go ??
		for ii=n:-1:1
            if nsteps[ii]>1
                maskf=mapslices(dvmaskexpand,maskf,dims = ii)
            end
        end
        maskc=maskf[coarsegridpoints...];
        ###### Test
        #maskc=trues(size(maskc))
        ######


        # Now scale pmn by the step factors

        pmnc=([ (1.0/nsteps[i])*pmn[i][coarsegridpoints...] for i=1:length(pmn) ]...,)

        # Check if Labs is a tuple of tuple; in this case also subsample

        if isa(Labs,Tuple)
            if isa(Labs[1],Number)
                Labsc=Labs;
            else
                Labsc=([ x[coarsegridpoints...] for x in Labs ]...,);
            end
        else
            # Create a tuple of L for the coarse grid; needed to be able to put some of them to zero
            #  Labsc=Labs;
            Labsc=(Labs*ones(n)...,);
        end


        # Create HI only of subsampled, otherwise use eye
        # Now prepare HI do go from the coarse grid to the fine grid. To do so
        # interprete de fine grid coordinates as those of pseudo-obs and use DIVAndtoos
        # Need the statevector strucure for the fine grid to go from grid to array

        svf = statevector_init((mask,))

        # For each coordinate in the tuplet xi, go from grid representation to tuplet
        # to have the pseudo-data coordinates


        Ic = localize_separable_grid(([statevector_pack(svf,(x,)) for x in xi]...,),maskc,xic);

        # Create fractional indexes of these data points in the coarse grid

        HI,outc,outbboxc = sparse_interp(maskc,Ic,iscyclic);
        HI = HI * sparse_pack(maskc)';

        ####################################
        # Advection constraint extracted
        ####################################



        otherargsc = Dict(otherargs)        

        # Search for velocity argument:
        if haskey(otherargsc,:velocity)
            # modify the parameter only in the coarse model
            vel = otherargsc[:velocity]
            otherargsc[:velocity] = ([ x[coarsegridpoints...] for x in vel]...,)
        end

        # For other constraints:
        # TODO TODO TODO TODO IF DIVAndjog should accept other constraints
        # Here should be straightfoward replace C by C*HI on the constraint structure





        # Prepare run of the coarse grid problem


        # For the coarse model, slightly adapth alphabc assuming a typical ratio of 4 is used
        # Search for alphabc argument:
        # modify the parameter only in the coarse model
        otherargsc[:alphabc] = 0.25

        # maybe try another norm using btrunc=3 or even btrunc=2 but full L here for the coarser version ?

        # Preconditionner with desactivated correlations in some directions

        # Hardcoded new test


		methodpccoarse=pcmethod

		tol=1e-3
		maxiter=20*Int(ceil(sqrt(size(HI)[1])))
		maxiter=minimum([2000,maxiter])
        diagshift=0.000001

		fi=0
		si=0

		#Method 1
		if methodpccoarse==1

			Labsccut=([Labsc[i]*lmask[i] for i=1:n]...,)
			# try classic 2D
			lmask1 = 0.0.*lmask;
			lmask1[1:2] .= 1.0;
			Labsccut=([Labsc[i]*lmask1[i] for i=1:n]...,)

			#fc,sc=DIVAndrun(maskc,pmnc,xic,x,f,Labsccut,epsilon2; otherargsc...,alpha=alphapc,btrunc=2)
			fc,sc=DIVAndrun(maskc,pmnc,xic,x,f,Labsccut,epsilon2; otherargsc...)

		    scP=1
			figuess=zeros(size(mask))
			xguess=1

			if !(sc==0)
                # TEST makes sure there are values in the coarse resolution solution
                # Take the fc coarse solution, pack to to statevector form
                # using sc.sv
                # Apply HI; this vector can also be used as a first guess for the PC

                # Preconditionner core
				scP=sc.P;
				xguess=HI*statevector_pack(sc.sv,(fc,));
				figuess,=statevector_unpack(svf,xguess)
			end
			# Try to clean up some memory here
			sc=0
			fc=0
			#gc()



			xr=randn(size(HI)[1],1)

			scalef=(xr'*(HI*(HI'*xr)))./(xr'*xr)
			scalef2=1 / scalef[1]

			scalef2=1.02
			xguess=xguess*scalef2
            figuess,=statevector_unpack(svf,xguess)

			# To compensate for the missing correlations in HI*scP*HI'

			diagshift=0.006*(sqrt(size(HI)[1]/size(HI)[2])-1);
			diagshift=0.03*(sqrt(size(HI)[1]/size(HI)[2])-1);
			diagshift=0.04*(sqrt(size(HI)[1]/size(HI)[2])-1);
			diagshift=0.02*(sqrt(size(HI)[1]/size(HI)[2])-1);

			#work1=0.*xguess # ::Array{Float64,1}
			#work2=0.*xguess # ::Array{Float64,1}
			#work3=0.*(HI'*xguess) # ::Array{Float64,1}
			work3=zeros(Float64,size(HI)[2])
			work3b=zeros(Float64,size(HI)[2])
			function compPC(iB,H,R)
			    #work1=zeros(Float64,size(HI)[1])
			    #work2=zeros(Float64,size(HI)[1])

			    #Possible improvement, instead of linear interpolation vi HI:
			    # HI' : unpack state vector, from the fine grid get grid points [coarsegridpoints...]
			    #  pack the coarse grid, apply scP, unpack to coarse grid, initialize a fine grid solution with this coarse grid solution
			    #  use DIVAnd_fill to fill the fine grid and finally pack the array
			    # maybe even accept a filtering filling ?
			    # But then with this approach the operator is not symmetric ?
                function fun!(x::Array{Float64,1},fx::Array{Float64,1})
					#@show size(x),typeof(HI),typeof(scP)
					#work1[:]=diagshift*x ::Array{Float64,1}

                    if VERSION >= v"0.7.0-beta.0"
					    mul!(work3::Array{Float64,1},(HI::SparseMatrixCSC{Float64,Int})',x::Array{Float64,1})
                    else
					    At_mul_B!(work3::Array{Float64,1},HI::SparseMatrixCSC{Float64,Int},x::Array{Float64,1})
                    end
					#work3[:]=HI'*x  ::Array{Float64,1}
					if scP==1
    					work3b=work3
					else
						work3b[:]=scP*work3 ::Array{Float64,1}
						# not yet defined for covariance
						#mul!(work3b::Array{Float64,1},scP,work3::Array{Float64,1})
					end

					#

					mul!(fx::Array{Float64,1},HI::SparseMatrixCSC{Float64,Int},work3b::Array{Float64,1})
					#fx[:]=HI*work3b ::Array{Float64,1}
					fx[:]=scalef2*fx ::Array{Float64,1}
					#fx[:]=work1+work2 ::Array{Float64,1}
					#mul!(fx,scalef2,fx)
					fx[:]=BLAS.axpy!(diagshift,x,fx)
                    #fx[:] = diagshift*x+scalef2*HI*(scP*(HI'*x));
                end

				#@time gc()
                return fun!
            end
			fi,si=DIVAndrun(mask,pmn,xi,x,f,Labs,epsilon2; otherargs...,
                            tol = tol,
                            maxit = maxiter,
                            inversion=:pcg,compPC = compPC, fi0 =figuess,btrunc=2)



		end

		#Method 2
		if methodpccoarse==2

			Labsccut=Labsc
			lmask1=0.0.*lmask
			lmask1[1:end] .= 1.0
			if n>3
				lmask1[3] = 0.0
			end
			Labsccut=([Labsc[i]*lmask1[i] for i=1:n]...,)
			fc,sc=DIVAndrun(maskc,pmnc,xic,x,f,Labsccut,epsilon2; otherargsc...,btrunc=3)
			scP=1
			figuess=zeros(size(mask))
			xguess=1
			if !(sc==0)
                # TEST makes sure there are values in the coarse resolution solution
                # Take the fc coarse solution, pack to to statevector form
                # using sc.sv
                # Apply HI; this vector can also be used as a first guess for the PC

                # Preconditionner core
				scP=sc.P;
				xguess=HI*statevector_pack(sc.sv,(fc,));
				figuess,=statevector_unpack(svf,xguess)
			end
			# Try to clean up some memory here
			sc=0
			fc=0
			#gc()

            xr=randn(size(HI)[1],1)

			scalef=(xr'*(HI*(HI'*xr)))./(xr'*xr)
			scalef2=1. / scalef[1]
			scalef2=1.
			xguess=xguess*scalef2

            figuess,=statevector_unpack(svf,xguess)


			# To compensate for the missing correlations in HI*scP*HI'


			diagshift=0.006*(sqrt(size(HI)[1]/size(HI)[2])-1);
			function compPC2cr(iB,H,R)
                function fun!(x,fx)
                    fx[:] = diagshift*x+scalef2*(HI*(scP*(HI'*x)));
                end
                return fun!
            end
			fi,si=DIVAndrun(mask,pmn,xi,x,f,Labs,epsilon2; otherargs...,
                            tol = tol,
                            maxit = maxiter,
                            inversion=:pcg,compPC = compPC2cr, fi0 =figuess,btrunc=2)


		end



		#Method 3
		if methodpccoarse==3
			lmask1=0.0.*lmask;
			lmask1[1:end] .= 1.0;
			if n>3
				lmask1[3]=0;
			end
			Labsccut=([Labsc[i]*lmask1[i] for i=1:n]...,)

			fc,sc=DIVAndrun(maskc,pmnc,xic,x,f,Labsccut,epsilon2; otherargsc...)
			PC2=1
			xguess=1
			if !(sc==0)
                # TEST makes sure there are values in the coarse resolution solution
                # Take the fc coarse solution, pack to to statevector form
                # using sc.sv
                # Apply HI; this vector can also be used as a first guess for the PC

                # Preconditionner core
				PC2=sc.P;
				xguess=statevector_pack(sc.sv,(fc,));
			end
			# Try to clean up some memory here
			sc=0
			fc=0
			#gc()

			PC1=1
			lmask1=0.0.*lmask;
			#lmask1[3:end]=1/1.42;
            #Try 3D
			if n>3
				lmask1[3]=1/1.42;
				Labsccut=([Labsc[i]*lmask1[i] for i=1:n]...,)

				fc,sc=DIVAndrun(maskc,pmnc,xic,x,f,Labsccut,10000.; otherargsc...)


				if !(sc==0)
					# TEST makes sure there are values in the coarse resolution solution
					# Take the fc coarse solution, pack to to statevector form
					# using sc.sv
					# Apply HI; this vector can also be used as a first guess for the PC

					# Preconditionner core
					PC1=sc.P;
				end
			    # Try to clean up some memory here
				sc=0
				fc=0
				#gc()
            end



			# tolerance on the gradient A x - b


			diagshift=0.001*(sqrt(size(HI)[1]/size(HI)[2])-1);


			xguess=PC1*(PC1*xguess);
			#xr=randn(size(HI)[1],1)
			#scalef=(xr'*(HI*(PC1*(PC2*(PC1*(HI'*xr))))))./(xr'*xr)
			xr=randn(size(HI)[2],1)
			scalef=(xr'*((PC1*(PC2*(PC1*(xr))))))./(xr'*xr)
			#xr=randn(size(HI)[2],1)
			scalefter=(xr'*(PC2*xr))./(xr'*xr)

			#Maybe adapt scaling including HI ?
			# xr=randn(size(HI)[1],1)
			# xz=HI'*xr;
			# scalef=(xz'*(PC1*(PC2*(PC1*xz))))./(xz'*xz)
			# xz=randn(size(HI)[2],1)
			# scalefter=(xz'*(PC2*xz))./(xz'*xz)


			scalef2=scalefter[1]/scalef[1]
			xguess=xguess*scalef2
			figuess,=statevector_unpack(svf,HI*xguess)

			function compPC3(iB,H,R)
                function fun!(x,fx)
                    fx[:] = diagshift*x+scalef2*(HI*(PC1*(PC2*(PC1*(HI'*x)))))
                end
                return fun!
			end

			fi,si=DIVAndrun(mask,pmn,xi,x,f,Labs,epsilon2; otherargs...,
                            tol = tol,
                            maxit = maxiter,
                            inversion=:pcg,compPC = compPC3, fi0 =figuess,btrunc=2)
		end



		if methodpccoarse==4

			Labsccut=Labsc
			lmask1=0.0.*lmask
			lmask1[1:end] .= 1.0;
			if n>3
				lmask1[3]=0.0;
			end
			Labsccut=([Labsc[i]*lmask1[i] for i=1:n]...,)
			fc,sc=DIVAndrun(maskc,pmnc,xic,x,f,Labsccut,epsilon2; otherargsc...,btrunc=3)
			scP=1


			PC2=1
			xguess=1
			if !(sc==0)
                # TEST makes sure there are values in the coarse resolution solution
                # Take the fc coarse solution, pack to to statevector form
                # using sc.sv
                # Apply HI; this vector can also be used as a first guess for the PC

                # Preconditionner core
				PC2 = sc.P
				xguess=statevector_pack(sc.sv,(fc,));
			end
			# Try to clean up some memory here
			sc=0
			fc=0
			#gc()



			#Now B2D on fine resolution model !
			Labsf=Labs
			if isa(Labs,Tuple)
                Labsf=Labs
            else
                Labsf=(Labs*ones(n)...,);
            end


			lmask1=0.0.*lmask;

			lmask1[1]=1.0/1.42;

			Labsccut=([Labsf[i]*lmask1[i] for i=1:n]...,)
			sc=0
			PC1a=1
			fc,sc=DIVAndrun(mask,pmn,xi,x,f,Labsccut,10000.; otherargs...,btrunc=3)
			if !(sc==0)
                # TEST makes sure there are values in the coarse resolution solution
                # Take the fc coarse solution, pack to to statevector form
                # using sc.sv
                # Apply HI; this vector can also be used as a first guess for the PC

                # Preconditionner core
				PC1a = sc.P
			end
			# Try to clean up some memory here
			sc=0
			fc=0
			#gc()

			lmask1=0.0.*lmask;

			lmask1[2]=1.0/1.42;
			Labsccut=([Labsf[i]*lmask1[i] for i=1:n]...,)
			sc=0
			PC1b=1
			fc,sc=DIVAndrun(mask,pmn,xi,x,f,Labsccut,10000.; otherargs...,btrunc=3)
			if !(sc==0)
                # TEST makes sure there are values in the coarse resolution solution
                # Take the fc coarse solution, pack to to statevector form
                # using sc.sv
                # Apply HI; this vector can also be used as a first guess for the PC

                # Preconditionner core
				PC1b = sc.P
			end
			# Try to clean up some memory here
			sc=0
			fc=0
			#gc()




			# tolerance on the gradient A x - b


			diagshift=0.006*(sqrt(size(HI)[1]/size(HI)[2])-1);


			xguess=PC1a*(PC1b*(PC1b*(PC1a*(HI*xguess))));
			xr=randn(size(HI)[1],1)
			scalef=(xr'*(PC1a*(PC1b*(HI*(PC2*(HI'*(PC1b*(PC1a*xr))))))))./(xr'*xr)
			xr=randn(size(HI)[2],1)
			scalefter=(xr'*(PC2*xr))./(xr'*xr)

			#Maybe adapt scaling including HI ?
			# xr=randn(size(HI)[1],1)
			# xz=HI'*xr;
			# scalef=(xz'*(PC1*(PC2*(PC1*xz))))./(xz'*xz)
			# xz=randn(size(HI)[2],1)
			# scalefter=(xz'*(PC2*xz))./(xz'*xz)


			scalef2=scalefter[1]/scalef[1]
			# Because of btrunc in major norm
			scalef2=scalef2*0.16*1.1
			xguess=xguess*scalef2
			figuess,=statevector_unpack(svf,xguess)


			function compPC4b(iB,H,R)
                function fun!(x,fx)
				    #return x -> diagshift*x-diagshift*(HI*(HI'*x))+scalef2*(HI*(PC1*(PC2*(PC1*(HI'*x)))))
				    #return x -> diagshift*x+scalef2*(PC1*(HI*(PC2*(HI'*(PC1*x)))))
				    #return x -> diagshift*x+scalef2*((HI*(PC2*(HI'*(x)))))
                    fx[:] = diagshift*x+scalef2*(PC1a*(PC1b*(HI*(PC2*(HI'*(PC1b*(PC1a*x)))))))
                end
                return fun!
			end

			fi,si=DIVAndrun(mask,pmn,xi,x,f,Labs,epsilon2; otherargs...,
                            tol = tol,
                            maxit = maxiter,
                            inversion=:pcg,compPC =compPC4b,fi0 =figuess,btrunc=2)


		end

	    if methodpccoarse==5
			lmask1=0.0.*lmask;
			lmask1[1:2] .= 1.0;
			Labsccut=([Labsc[i]*lmask1[i] for i=1:n]...,)

			fc,sc=DIVAndrun(maskc,pmnc,xic,x,f,Labsccut,epsilon2; otherargsc...)
			PC2=1
			xguess=1
			if !(sc==0)
                # TEST makes sure there are values in the coarse resolution solution
                # Take the fc coarse solution, pack to to statevector form
                # using sc.sv
                # Apply HI; this vector can also be used as a first guess for the PC

                # Preconditionner core
				PC2=sc.P;
				xguess=statevector_pack(sc.sv,(fc,));
			end
			# Try to clean up some memory here
			sc=0
			fc=0
			#gc()
			lmask1=0.0.*lmask;
            lmask1[3:end] .= 1.0/1.42
			Labsccut=([Labsc[i]*lmask1[i] for i=1:n]...,)
			PC1=1
			sc=0

			if size(lmask)[2]> 2
				fc,sc=DIVAndrun(maskc,pmnc,xic,x,f,Labsccut,10000.; otherargsc...)
			end

			if !(sc==0)
                # TEST makes sure there are values in the coarse resolution solution
                # Take the fc coarse solution, pack to to statevector form
                # using sc.sv
                # Apply HI; this vector can also be used as a first guess for the PC

                # Preconditionner core
				PC1=sc.P;
			end
			# Try to clean up some memory here
			sc=0
			fc=0
			#gc()

            #Now B2D on fine resolution model !
			if isa(Labs,Tuple)
                Labsf=Labs
            else
                Labsf=(Labs*ones(n)...,);
            end
			lmask1=0.0.*lmask;
			lmask1[1]=1.0/1.42;
			Labsccut=([Labsf[i]*lmask1[i] for i=1:n]...,)
			sc=0
			PC1a=1
			fc,sc=DIVAndrun(mask,pmn,xi,x,f,Labsccut,10000.; otherargs...)
			if !(sc==0)
                # TEST makes sure there are values in the coarse resolution solution
                # Take the fc coarse solution, pack to to statevector form
                # using sc.sv
                # Apply HI; this vector can also be used as a first guess for the PC

                # Preconditionner core
				PC1a = sc.P
			end
			# Try to clean up some memory here
			sc=0
			fc=0
			#gc()

			lmask1=0.0.*lmask;
			lmask1[2]=1.0/1.42;
			Labsccut=([Labsf[i]*lmask1[i] for i=1:n]...,)
			sc=0
			PC1b=1
			fc,sc=DIVAndrun(mask,pmn,xi,x,f,Labsccut,10000.; otherargs...)
			if !(sc==0)
                # TEST makes sure there are values in the coarse resolution solution
                # Take the fc coarse solution, pack to to statevector form
                # using sc.sv
                # Apply HI; this vector can also be used as a first guess for the PC

                # Preconditionner core
				PC1b = sc.P
			end
			# Try to clean up some memory here
			sc=0
			fc=0
			#gc()



			# tolerance on the gradient A x - b


			diagshift=0.005*(sqrt(size(HI)[1]/size(HI)[2])-1);
			diagshift=0.00015*(sqrt(size(HI)[1]/size(HI)[2])-1);

			xguess=(PC1b*(PC1a*(PC1a*(PC1b*(HI*(PC1*(PC1*xguess)))))));
			xr=randn(size(HI)[1],1)
			scalef=(xr'*(PC1a*(PC1b*(HI*(PC1*(PC2*(PC1*(HI'*(PC1b*(PC1a*xr))))))))))./(xr'*xr)
			#xr=randn(size(HI)[2],1)
			#scalef=(xr'*((PC1*(PC2*(PC1*(xr))))))./(xr'*xr)
			xr=randn(size(HI)[2],1)
			scalefter=(xr'*(PC2*xr))./(xr'*xr)

			#Maybe adapt scaling including HI ?
			# xr=randn(size(HI)[1],1)
			# xz=HI'*xr;
			# scalef=(xz'*(PC1*(PC2*(PC1*xz))))./(xz'*xz)
			# xz=randn(size(HI)[2],1)
			# scalefter=(xz'*(PC2*xz))./(xz'*xz)


			scalef2=scalefter[1]/scalef[1]


			scalef2=scalef2
			xguess=xguess*scalef2
			figuess,=statevector_unpack(svf,xguess)

			function compPC5(iB,H,R)
                function fun!(x,fx)
				    #return x -> diagshift*x-diagshift*(HI*(HI'*x))+scalef2*(HI*(PC1*(PC2*(PC1*(HI'*x)))))
				    fx[:] = diagshift*x+scalef2*(PC1a*(PC1b*(HI*(PC1*(PC2*(PC1*(HI'*(PC1b*(PC1a*x)))))))))
                end
                return fun!

			end

			fi,si=DIVAndrun(mask,pmn,xi,x,f,Labs,epsilon2; otherargs...,
                            tol = tol,
                            maxit = maxiter,
                            inversion=:pcg,compPC = compPC5, fi0 =figuess,btrunc=2)


            # posteriori=fi.*figuess
			# Some ideas for error calculations
			#fi,si=DIVAndrun(mask,pmn,xi,x,f,Labs,epsilon2; otherargs...,pcargs...,inversion=:pcg,compPC = DIVAnd_pc_sqrtiB, fi0 =xguess)
			#errfield=diagMtCM(sc.P,HI')
			#erri,=statevector_unpack(si.sv,errfield)
			# For error field based on coarse one, use DIVAnd_filter3 with ntimes=Int(ceil(mean(nsteps)))
			# First test
			# Possible optimization for a climatology production, for each tile return diagshift and niter with some random changes in diagshift and search
			# for optimum. In particular if several climatology runs are produced (slight changed parameters), the preliminary runs, without error calculations for example
			# could provide good estimates ?

		end

		# A posteriori scaling of initial guess
		#	fs=statevector_pack(si.sv,(fi,));
		#	fgs=statevector_pack(si.sv,(figuess,));
		#	al=dot(fs,fgs)/dot(fgs,fgs)
		#	@show al,dot(fgs-fs,fgs-fs)/dot(fs,fs),dot(al*fgs-fs,al*fgs-fs)/dot(fs,fs)



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

# LocalWords:  fi DIVAnd pmn len diag CovarParam vel ceil moddim fracdim
