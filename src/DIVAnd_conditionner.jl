"""

    Lpmnmean = DIVAnd_Lpmnmean(pmn,len);

# In each direction, calculates the mean value of the length scale times the metric in this direction
# So it basically looks for the average resolution on the grid

# Input:

* `pmn`: scale factor of the grid. pmn is a tuple with n elements. Every
       element represents the scale factor of the corresponding dimension. Its
       inverse is the local resolution of the grid in a particular dimension.

* `len`: correlation length


# Output:

* `Lpmnmean`: Array of mean values of L times metric

"""
function DIVAnd_Lpmnmean(pmn::NTuple{N,Array{T,N}}, len) where {N,T}
    Lpmnmean = Vector{Float64}(undef, N)

    for i = 1:N
        if isa(len, Number)
            Lpmnmean[i] = mean(len * pmn[i])
        elseif isa(len, Tuple)

            if isa(len[1], Number)
                Lpmnmean[i] = mean(len[i] * pmn[i])
            else
                Lpmnmean[i] = mean(len[i] .* pmn[i])
            end

        end
    end

    return Lpmnmean
end


function OIcorrregijk!(corr,Lpmnmean;kernelfun=(x-> (1 .+x).*exp.(-x)))   #(x<=0 ? 1.0 : sqrt.(x).*besselk.(1,sqrt.(x)))))
    # returns array corr of correlations assuming a regular grid and uniform L
    # Will speed up subsequent OI interpolation
    # kernelfunction in r where r is the nondimensional distance.
    # 
    R=CartesianIndices(corr)
    I1=oneunit(first(R))
    
    # 
    LL=[0.50765164, 1.0, 1.370937, 1.677412167, 1.943162456]
    
    n=ndims(corr)
    m=2
    if n>2
        m=3
    end
    twonu=2*m-n
    
    if mod(twonu,2)==0
        kernelfun=x->(x<=0 ? 1.0 : x.*besselk.(1,x))
    end
    iLS=LL[Int(twonu)]./Lpmnmean
    for I in R
        corr[I]=0.0
        dists=0.0
        for ij=1:ndims(corr)
            dx=(I[ij]-I1[ij])*iLS[ij]
            dists=dists+dx*dx
        end
        corr[I]=kernelfun(sqrt(dists))
        
    end
    
end


# cholesky decomposition of HBH'+R 
function DIVAndOIBRdecomposition(corr,posdata,epsilon2)
    BPR=ones(Float64,size(posdata,1),size(posdata,1))
    I1=oneunit(posdata[1])
    
    for i=1:size(posdata,1)
        for j=i+1:size(posdata,1)
            work=posdata[i]-posdata[j]
            work=max(work,-work)+I1
            BPR[i,j]=corr[work]
            BPR[j,i]=BPR[i,j]
        end
        BPR[i,i]=BPR[i,i]+epsilon2[i]
    end
    return factorize(Symmetric(BPR))
end



# Diagonal preconditionning
function compPCDIAG(iB,H,R)
    hdiag=zeros(Float64,size(iB,1))
    for i=1:size(iB,1)
        hdiag[i]=1.0/(iB[i,i]+H[:,i]'*(R\H[:,i]))
    end
    # Warning, if btrunc is used, iB is not complete
    function fun!(x,fx)
        
        
            
            fx.=x.*hdiag
        
       
    end
    return fun!

    
end


#### TO CLEAN UP AND OTPIMIZE


function DIVAndFFTpcg(mask,pmn,xyi,xy,f,len,epsilon2;moddim=zeros(Int32,ndims(mask)),maxsuper=-1,onlyB=false)

# For superobs: find closest cartesianindex (using I = localize_separable_grid   and make list.
# Check that list has no duplicates (could arrive if a lot of points are clustered)
# Create cholesky(HBH'+R)  and story it for further use
# Initial guess from OI
# Then create precontioning function fun!


    
    corr=zeros(Float64,size(pmn[1]))
    Lpmnmean=DIVAnd_Lpmnmean(pmn,len)
    OIcorrregijk!(corr,Lpmnmean)
    
# ideally 1.2 superobs in each influence bubble domain. Maybe adapt in higher dimensions (sphere vs cube)
# roughly factor 2 in 3D
# and look at spread of data points?
# so (max(x)-min(x))/mean(lenx)
    super=maxsuper
    if super<0
    super=2*Int(ceil(1.2^ndims(pmn[1])*prod(size(pmn[1])./Lpmnmean)))
    end
	super=min(super,size(f,1))
    @show super
# TODO, use proper weighting if epsilon2 is not constant
    newx,newval,sumw,varp,idx=DIVAnd_superobs(xy,f,super)

    ILOC = localize_separable_grid(newx,mask,xyi)
      myval=[]
      myloc=[]
      myeps=[]
      for ii=1:size(ILOC,2)
      if sum(ILOC[:,ii].>0)==size(ILOC,1)
         myloc=[myloc...,CartesianIndex(Int.(round.(ILOC[:,ii]))...)]
         myval=[myval...,newval[ii]]
         myeps=[myeps...,1/sumw[ii]]
         else
        
      end
      end
      if size(unique(myloc))==size(myloc)
        @show "seems ok",size(myloc)
          else
        # To do: regroup points that are collocated.
        @show "NEED TO ADD DEALING WITH SUPERPOSED SUPEROBS"
        # Sort superobs by Indices
        # myloclinear=LinearIndices(mask)[[myloc...]]
        # Then go through the list and drop points.
        # 
        # Then back to Cartesian
        # myloc=CartesianIndices(mask)[[myloclinear...]]
      end
      #@show epsilon2,myeps,myloc
      BRD=DIVAndOIBRdecomposition(corr,myloc,epsilon2.*myeps)
      xOI= BRD\myval
      fOI=deepcopy(xyi[1])
    
    
    
    

    
    
@show (2 .* ones(Int32,ndims(corr)).- moddim) .* size(corr)
gaussf=zeros(Float64,((2 .* ones(Int32,ndims(corr)) .- moddim) .*size(corr)...))



gaussf[CartesianIndices(corr)].=corr

function halfsize(x)
    if mod(x,2)==0
        return Int(x/2)
    else
        return Int((x+1)/2)
    end
end

for i=1:ndims(gaussf)
    @show i
    idxt = ntuple(x->x == i ? (size(gaussf,i):-1:size(gaussf,i)-halfsize(size(gaussf,i))+2 ) : Colon(), ndims(gaussf))
    idxf = ntuple(x->x == i ? (2:halfsize(size(gaussf,i)) ) : Colon(), ndims(gaussf))
    @show idxf,idxt,size(gaussf)
    gaussf[idxt...]=gaussf[idxf...]
end

    PFFT=plan_rfft(gaussf)
    FFTG=PFFT*gaussf

    PiFFT=plan_irfft(FFTG,size(gaussf,1))

    xval=zeros(Float64,size(gaussf))
    
    work=zeros(Float64,size(myloc,1))
    
    #gaussf=[]
	# reuse corr and gaussf
    xin=corr
    ork=gaussf #zeros(Float64,size(gaussf))
	# probably not needed but anyway
    GC.gc()


    
        xin[myloc].=xOI
        xval[CartesianIndices(xin)].=xin
        ork .=PiFFT*((PFFT*xval).*FFTG)
        fOI.=ork[CartesianIndices(xin)]
    

      # Maybe increase value when grid is stretched ?
      mymax= Int(round(1.2*sqrt(prod(size(mask)))))
	  # If only preconditionned by B, increase
      if onlyB
	  mymax=mymax*1.5
	  end
     #
  function compPCFFT(iB,H,R)
  
    function fun!(x,fx)
        
        #  avoid allocation
        xin[:],=statevector_unpack(mysv, x)
        xval[CartesianIndices(xin)].=xin
        ork .= PiFFT*((PFFT*xval).*FFTG)
        fx[:].=statevector_pack(mysv, (ork[CartesianIndices(xin)],))
              
        # if return just Bx
        if onlyB
     		return
		end
		# ork is Bx 
        work[:].=BRD\ork[myloc]
        xin[:].=0.0
        xin[myloc].=work
        xval[CartesianIndices(xin)].=xin
        ork .=PiFFT*((PFFT*xval).*FFTG)
        fx[:].= fx.-statevector_pack(mysv, (ork[CartesianIndices(xin)],))
        
    end
    return fun!

  end
    
	
  return fOI,compPCFFT,mymax
end





# Copyright (C) 2008-2023 Alexander Barth <barth.alexander@gmail.com>
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
