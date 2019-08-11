"""
Computes a  heatmap based on locations of observations

dens2 = DIVAnd_heatmap(mask,pmn,xi,x,inflation,Labs;Ladaptiveiterations=0,myheatmapmethod="DataKernel",
    optimizeheat=true,otherargs...)


# Input:
*  `mask`: mask as usual
*  `pmn` : tuple of metrics as usual
*  `xi`: tuple of coordinates of the grid for the heatmap
*  `x` : tuple of coordinates of observations
*  `inflation`: array generally of ones. For some applications an observation can carry a different weight which is then encoded in the as
*  `Labs` : the length scales for diva. Here their meaning is the spread of the observations for the Kernel calculation
*              if zero is provided, the routine applies an empirical estimate

*   `Ladaptiveiterations`: adaptive scaling where the length scales are adapted on the data density already estimated. You can iterate
*   `optimizeheat` : boolean which can turn on or off an algorithmic optimisation. Results should be identical
*   `myheatmapmethod`: can be "Automatic", "GridKernel" or "DataKernel". 

*   `otherargs...`: all other optional arguments DIVAndrun can take (advection etc)

# Output:
*  `lambda`: data density field (integral is one)
"""
function DIVAnd_heatmap(mask,pmn,xi,x,inflation,Labs;Ladaptiveiterations=0,myheatmapmethod="DataKernel",
    optimizeheat=true,otherargs...)
#
    
    # Create output array on the same grid as mask pmn and xi
    dens2=zeros(Float64,size(mask))
    #dens2opt=zeros(Float64,size(mask))
    # Dimensionality of the problem
    DIMS=ndims(mask)
    NP=size(inflation)[1]
#Empirial estimate Silverman's (1986) rule of thumb
    LHEAT=Labs
    
    mymethod=myheatmapmethod
	
	if myheatmapmethod=="Automatic"
		mymethod="DataKernel"
		
		if NP>sum(mask.==true)
		mymethod="GridKernel"
		end
	
	
	end
	
	
    trytooptimize=optimizeheat
    #mymethod="GridKernel"
# 
    if Labs==0
        # Estimate
    
    
        varx=zeros(Float64,DIMS)
        LF=zeros(Float64,DIMS)
        
        for i=1:DIMS
            meanxo=sum(inflation.*x[i])/sum(inflation)
            varx[i]=sum(inflation.*(x[i].-meanxo).^2)/sum(inflation)
            #@show i,meanxo,varx[i]
            LF[i]=sqrt(varx[i])/((DIMS+2.0)*NP/4.0)^(1.0/(4.0+DIMS))
        end
    

        LHEAT=(LF...,)
        @show "Estimation", LHEAT
        
    end
    
    
    #  Option: Make iterations to include adaptive L where L is readjusted based on the calculated density
    #  Contrary to classical KDE no need to chose between balloon or the other method as L is now in the correlation, not in the scaling of the 
    # "distance". Probably much better
    
    # 
    Lfortuple=Array{Any}(undef,DIMS)
    Ltuple=LHEAT
    
    #Ltuple=(Lfortuple...,)
    
#     
# Co
# 
    inflationsum=0
    
    
   
    
    
    for Literations=1:1+Ladaptiveiterations
        
        inflationsum=0
        dens2 .= 0.0
        #dens2opt .= 0.0
        
        xxx=Array{Any}(undef,DIMS)
        if trytooptimize
        #@show "Try to calculate a decomposition" 
            #Decompose once and for all
			if mymethod=="DataKernel" 
            FIopt,Sopt=DIVAnd.DIVAndrun(mask,pmn,xi,x,inflation,Ltuple, 1.0E10 ;otherargs...)
            svf=statevector_init((mask,))
            #@show "a",size(Sopt.H'),
            #eiarr=zeros(size(inflation))
            #@time vv= Sopt.P.factors.PtL \ (Sopt.H'*eiarr)
            #@show "b"
            #@time vb=Sopt.P.factors.UP \ vv
            
            #@show size(vb)
            #ongrid,=statevector_unpack(svf,vb)
            #@show size(ongrid)
            #@show vv*inflation,size(inflation)
			end
			if mymethod=="GridKernel"
			 FIopt,Sopt=DIVAnd.DIVAndrun(mask,pmn,xi,x,inflation,Ltuple, 1.0E10 ;otherargs...)
            svf=statevector_init((mask,))
			end
        end
        
        # VERSION A: covariance of one data point with grid points
        if mymethod=="DataKernel"
            
        for myi=1:NP
        
            if trytooptimize
                eiarr=zeros(size(inflation))
                eiarr[myi]=1
                vv= Sopt.P.factors.PtL \ (Sopt.H'*eiarr)
                vb=Sopt.P.factors.UP \ vv
                   
                FI,=statevector_unpack(svf,vb)
                integ=DIVAnd_integral(mask,pmn,FI)
                
                    
            else    
                
                
            for ii=1:DIMS
                xxx[ii]=[x[ii][myi]]
            end
            #@show xxx,LF,size(xi[1]),size(pmn[1])
              #  Use of WOODBURY and decomposed B in s could make it faster
                FI,S=DIVAnd.DIVAndrun(mask,pmn,xi,(xxx...,),[1.0],Ltuple,0.001;otherargs...)
        
    # Add here the constraint that each integral is one: accepts non unit values vi inflation to reflect mulitple observations
    # Also accept errors on obs? But how ? In the integral so it has less influence on overall sum (which will be scaled again?)
    #
                integ=DIVAnd_integral(mask,pmn,FI)
        #@show integ
                    
        end            
            if integ !=0
                dens2=dens2 .+ inflation[myi]*FI/integ
                inflationsum=inflationsum+inflation[myi]
            else
                @show "?? Not on grid ?",integ,sum(FI),xxx
            end
           
           
                
        end
        
        
            dens2=dens2/inflationsum
            
            #@show sum(dens2),DIVAnd_integral(mask,pmn,dens2),NP,inflationsum
            
        end
        # VERSION B: covariance of one grid point with all data points
        # TODO: implement optimized version .....
        if mymethod=="GridKernel"
		
		  if trytooptimize
		    svsize=sum(mask.==true)
			xdens=zeros(Float64,svsize)
			xval=zeros(Float64,svsize)
		    @show "OPti",svsize
			for myi=1:svsize
			eiarr=zeros(Float64,svsize)
            eiarr[myi]=1.0
			vv= Sopt.P.factors.PtL \ eiarr
            vb=Sopt.P.factors.UP \ vv
                
                FI,=statevector_unpack(svf,vb)
                integ=DIVAnd_integral(mask,pmn,FI)
				#@show integ
			vb=vb/integ
			#@show size(vb),size(inflation),size((Sopt.H*vb))
			xdens[myi]=sum(inflation .* (Sopt.H*vb))
				
			end
			dens2,=statevector_unpack(svf,xdens)
		    dens2=dens2/DIVAnd_integral(mask,pmn,dens2)
		  
		  else
            xaugmented=Array{Any}(undef,DIMS)
            @show "Grid based Kernels" 
            Rinf=deepcopy(inflation)
            Rinf.=1.0E10
            Raugmented=[Rinf...,0.000001]
            valaugmented=[inflation...,1.0]
            for myi in eachindex(mask)
                
                if mask[myi]
                    # Only on grid:
                    # Add one virtual data point which is the grid point. 
                    for ii=1:DIMS
                        xaugmented[ii]=[x[ii]...,xi[ii][myi]]
                    end
                    
                    # The real data points with infinite R
                    # Use of WOODBURY and decomposed B in s could make it faster
                    FI,s=DIVAnd.DIVAndrun(mask,pmn,xi,(xaugmented...,),valaugmented,Ltuple,Raugmented;otherargs...)
                    
                    # Calculate normalization constant
                    integ=DIVAnd_integral(mask,pmn,FI)
                    # After analysis, retrieve via S the analysis at those points and apply normalization constant
                    # as well as inflation
                    bidon=(s.obsconstrain.H)*statevector_pack(s.sv,(FI,))
                    dens2[myi]=sum(bidon[1:end-1].*inflation)/integ
                    
                    #@show myi,integ
                
                    
                end
            
            end
            # Need for overall scaling ?
            dens2=dens2/DIVAnd_integral(mask,pmn,dens2)
          end
		end
        
        
        
        
        
        
        
        
        # Optional L scaling

        if Ladaptiveiterations>0
            @show Literations,size(Ltuple[1])
            lambda=DIVAnd_scaleL(mask,pmn,dens2)
            for jj=1:DIMS
                @show LHEAT[jj]
                Lfortuple[jj] = LHEAT[jj] .* lambda
            end
            Ltuple=(Lfortuple...,)
        end
        
        
        

    end    
        
        
    return dens2
    
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


