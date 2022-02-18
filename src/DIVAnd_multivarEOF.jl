function DIVAnd_multivarEOF(mask,pmn,xi,x,f,lenin,epsilon2in;eof=(),velocity=(),kwargs...)


"""
    DIVAnd_multivarEOF(mask,pmn,xi,x,f,len,epsilon2; <keyword arguments>)

Perform an n-dimensional multivariate variational analysis of the observations `f` located at
the coordinates `x`. The output array `fi` represent the interpolated field at the grid
defined by the coordinates `xi` and the scales factors `pmn`.
The dimensions of the variables are arbitrary (generally space and time coordinates), but the LAST dimension must be the dimension for the different variables.
So for a multivariate analysis with 3 variables, the last dimension has three components and the corresponding coordinates are just 1, 2, 3 (the index of the variable).
All coordinates, including the input data needs to apply the same rule, so if your normal data have (x,y,z) coordinates for example, you will now have (x,y,z,v) where v contains for each data point the variable
to which it corresponds (1, 2, 3 ...)


# Input: same as DIVAndrun
* `mask`: binary mask delimiting the domain. true is inside and false outside.
For oceanographic application, this is the land-sea mask where sea is true and land is false.

* `pmn`: scale factor of the grid. pmn is a tuple with n elements. Every
   element represents the scale factor of the corresponding dimension. Its
   inverse is the local resolution of the grid in a particular dimension.
   For example, in two dimensions, `pmn` is a tuple `(pm,pn)` where `pm` is
   the inverse of the local resolution in first dimension and `pn` is the the inverse
   of the local resolution in second dimension.

*  `xi`: tuple with n elements. Every element represents a coordinate
  of the final grid on which the observations are interpolated.

* `x`: tuple with n elements. Every element represents a coordinate of
  the observations.

* `f`: value of the observations *minus* the background estimate (vector of
  `m` elements where `m` is the number of observations). See also note.

* `len`: tuple with n elements. Every element represents the correlation length for a given dimension.
     HERE THE LAST DIMENSION MUST HAVE ZERO L

* `epsilon2`: error variance of the observations (normalized by the error variance of the background field). `epsilon2` can be a scalar (all observations have the same error variance and their errors are decorrelated), a vector (all observations can have a different error variance and their errors are decorrelated) or a matrix (all observations can have a different error variance and their errors can be correlated). If `epsilon2` is a scalar, it is thus the *inverse of the signal-to-noise ratio*.

# Optional input arguments specified as keyword arguments and can take the same arguments as DIVAndrun.     
#    One optional keyword argument

*  `eof`: if provided contains an array of the
       coefficients from a linear fit between the variables (typically obtained by
       preliminary statistics or a run of this multivariate routine on a larger data set
       and which produced the eof coefficients). `eof`=[1.0,1.0] for example means the two variables are positively correlated with a slope 1.
	   `eof`=[1.0,-0.5] means the variables are negatively correlated and that for a variable 1 value of 1, variable 2 is expected to have a value of -0.5



#  Output: 

* `fi`:  multivariate field with the last dimension corresponding to the different variables
* `s`:    structure as for divandrun
* `eof`:  eof between variables (see explanation for the input keyword)
* `eofamplitudes`:   analysed amplitudes of the eofs. If "multiplied" by the eof that provides a new background field, a "common field between variables" which was subtracted for the classical analysis
* `emap`:  error map from univariate approach
* `emapm`: error map from multivariate approach 

"""

    
    function reducedims(xx)
    # utility function to prepare analysis on a grid WITHOUT the last dimension
    # so need to drop last element of the tuple xx and in each remaining element of the tuple the last dimension. 
    # Here done by summing (in principle same values but who knows how this will be used one day)    
          return tuple([dropdims(sum(xx[j],dims=ndims(xx[1]))/size(xx[1])[end],dims=ndims(xx[1])) for j=1:ndims(xx[1])-1]...)
    end
    
    # Later we need data individual epsilon2, so we enforce them here to start with
    nd=size(x[1])[1]
    epsilon2=deepcopy(epsilon2in)
    if isa(epsilon2in,Number)
        epsilon2=epsilon2in.* ones(nd)
    end
    
       
    
    # Some saveguards here; len must be zero on last dimension
    # At the same time create the lenlow correlation length for the analysis
    # in the lower dimension (without "variable" dimension) for the eof amplitudes
    #
    len=deepcopy(lenin)
    lenlow=[]
    eofin=eof
   
        if isa(lenin, Number)
            @warn "Expecting  len array with zero in last direction; will be enforced"
            lenlow=tuple(repeat([lenin],ndims(mask)-1)...)
            len=tuple(repeat([lenin],ndims(mask)-1)...,0.0)
            @debug "Len enforced $len"
        elseif isa(lenin, Tuple)

            if isa(lenin[1], Number)
                lenlow=(len[1:end-1])
            
               if lenin[end]!=0.0
                @warn "Expecting zero len in last direction; enforcing" 
                #@show lenin[end]
                #@show lenin[1:end-1]
                len=tuple(lenin[1:end-1]...,0.0)
                @debug "len enforced, $len"
                end
            else
                lenlow=reducedims(lenin)
                if sum(len[end])!=0.0
                @warn "Expecting zero len in last direction , enforcing"   
                len=(len[1:end-1]...,0.0 .* len[1])
                end
            end

        end
  
    
    
    # also need to deal with advection constraint dimension reduction
    velocitylow=()
    if velocity != ()
        velocitylow=reducedims(velocity)
    end
    # Dimension reduction
    masklow=reshape(mask[1:prod(size(mask)[1:end-1])],size(mask)[1:end-1])
    pmnlow=reducedims(pmn)
    xilow=reducedims(xi)
    xlow=(x[1:end-1])
    orivar=var(f)
    # Variance per layer (i.e. per variable)
    orivarl=zeros(Float64,size(mask)[end])
    for ll=1:size(mask)[end]
        orivarl[ll]=var(f[abs.(x[end].-ll).<0.3])
    end
    
    
    # univariate analysis 
    fi,s=DIVAndrun(mask,pmn,xi,x,f,len,epsilon2;velocity=velocity,kwargs...)        
    emap,meth=DIVAnd_errormap(mask,pmn,xi,x,f,len,epsilon2,s;method="cheap",velocity=velocity,kwargs...)        
    @debug "error method in multivar $meth"
    # keep points where error is small enough
    limite=0.5
    whereeof=dropdims(sum(emap,dims=ndims(fi)).<limite*size(mask)[end],dims=ndims(fi))
    
    eofamplitudes=zeros(Float64,size(mask)[1:end-1])
    amplitudes=zeros(Float64,size(f))
    work=zeros(Float64,size(fi)[ndims(fi)],sum(whereeof))
    if sum(whereeof)<30
               @warn "Sorry not reliable to infer correlations"
               return fi,s,eof,eofamplitudes,emap,emap
    end
    
    # add a safeguard based on number of observations in each layer
    # but only if eof calculation is demanded
    if eofin==()
      notenough=false
      
      for ll=1:size(mask)[end]
        if sum(abs.(x[end].-ll) .< 0.5 ) < 30
            notenough=true
            @warn "not enough data in multivariate layer $ll"
        end
      end
      if notenough
        return fi,s,eof,eofamplitudes,emap,emap
      end
    end
    # End check on amount of data
    
    
        # Now loop on improving the univariate approach
        # fixed number for the moment
        ML=1
        if eofin==()
           ML=3
        end
        for mloop=1:ML
            
            # Collect locations where errors are low in all fields 
            # to allow for eof calculations
             lambda=[]
             if eofin==()            
               for i=1:size(fi)[ndims(fi)]
                work[i,:]=selectdim(fi, ndims(fi), i)[whereeof]
               end
             
            # Now calculate EOF
              U,lambda,V=svd(work)
              @debug lambda
              eof=U[1,:]
              @debug eof
            
              if  lambda[1]-lambda[2] < lambda[2]/2
               @warn "Sorry too close eigenvalues $lambda"
                      return fi,s,eof,eofamplitudes,emap,emap
              end
            
                   else
               # eof was provided as input
               # so just scale it 
               eof.= eof / norm(eof)
               # and provide estimate of eigenvalues ?
               # TODO try to find a proxy
               lambda=ones(Float64,size(eof)[1])
             end
        # Check for eof not including weak component
        
        for k=1:size(eof)[1]
                if abs(eof[k])<0.1/size(eof)[1]
                      @warn "Sorry too weak  correlations layer $k"
                      return fi,s,eof,eofamplitudes,emap
                end
        end
        
        
            
            
            # now project data to get amplitudes
            for j=1:size(amplitudes)[1]
                k=round(Int,x[ndims(fi)][j])
                amplitudes[j]=f[j]/eof[k]#/(1+epsilon2)
            # regularize expliting lamdas ?
            end
            @debug extrema(amplitudes)
            # now make analysis of the amplitudes with one dimension less
            # and take into account that you only have part of the signal in the
            # amplitudes. Ad hoc formula for the moment.
            # 
            
            epsilon2low=epsilon2 .+ sum(lambda[2:end].^2)/lambda[1]^2
            eofamplitudes,sa=DIVAndrun(masklow,pmnlow,xilow,xlow,amplitudes,lenlow,epsilon2low;velocity=velocitylow,kwargs...)
            
            
            @debug size(eofamplitudes)
            @debug size(eof[:,:])
            # from amplitudes and eof construct full field
            eofbasedfield=reshape(reshape(eofamplitudes,:)*eof[:,:]',size(mask))
            
            # now amplitudes will keep residuals
            
            amplitudes .= s.obsconstrain.yo-
           (s.obsconstrain.H)*statevector_pack(s.sv, (eofbasedfield,))
            
                   
            @debug orivar, var(amplitudes)
         
            # For layers with less information left increase more 
            for ll=1:size(mask)[end]
                    wh=abs.(x[end].-ll).<0.3
                    epsilon2low[wh].=epsilon2[wh].*orivarl[ll]./var(amplitudes[wh])
            end
        
            # Now try to squeeze out information from the residuals (not represented by eof)
            ww,sa=DIVAndrun(mask,pmn,xi,x,amplitudes,len,epsilon2low;velocity=velocity,kwargs...)
            # now add indidivual contributions for full analys
            fi=eofbasedfield+ww
            
                 
        
        end 
    # Once the iterations finished
    #For the error map, include correlation by adding virtual points in each layer ?
    # Add hoc version:
    # each point of another layer has provided an information, but not as much as a point from the layer itself
    # the ratio of information is probably ... hence pseudo epsilon2 =
      nlay=size(mask)[end]
      nd=size(x[1])[1]
    # Arrays to store pseudo data
      xm=[[repeat(x[j],nlay) for j=1:ndims(mask)]...]
      
    # Do this also by layer values of ?
    # maybe fill in with inflated values see before then put
    # original values only for original points back
    # TODO
      #inflation=orivar/var(amplitudes)
      #@debug inflation,size(xm),size(epsilon2m)
      # epsilon for a points used in other layers
       epsilon2mm=deepcopy(epsilon2)
       for ll=1:size(mask)[end]
        www=abs.(x[end].-ll).<0.3
        epsilon2mm[www].=epsilon2mm[www] .*orivarl[ll]./var(amplitudes[www])
       end
      # copy for each layer
      epsilon2m=repeat([epsilon2mm...],nlay)
      # but replace by original  for real layer data
      epsilon2m[1:size(epsilon2)[1]].=epsilon2
      for k=2:nlay
        coo=ndims(mask)
        # just take the points from one layer and copy them into other layers
		# add 0.0001 to make sure rounding is not a problem
        xm[coo][(k-1)*nd+1:k*nd].= mod.(xm[coo][(k-2)*nd+1:(k-1)*nd] .+0.0001,nlay).+1
      end
      # Values of data are unimportant for the error field. So just repeated
      fm=repeat(f,nlay)
    
      # Error maps for the multivariate approach
      emapm,methm=DIVAnd_errormap(mask,pmn,xi,tuple(xm...),fm,len,epsilon2m,s;method="cheap",velocity=velocity,kwargs...)  
      @debug methm
    # Banzaii, finished      
    return fi,s,eof,eofamplitudes,emap,emapm
end