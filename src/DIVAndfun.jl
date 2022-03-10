

function DIVAndfun(x,f;mask=nothing,pmn=nothing,xi=nothing,len=nothing,epsilon2=nothing,kwargs...)
    
    
    
"""
	
	myfunny=DIVAndfun(x,f;mask=nothing,pmn=nothing,xi=nothing,len=nothing,epsilon2=nothing,kwargs...)
	
Provides a simplified interface to  DIVAndrun and in return provides an interpolation FUNCTION rather than the gridded field of DIVAndrun

You can use ALL paramters you can use with DIVAndrun, but some of them are made optional here by using keywords.

The necessary input is the tuple of coordinates at which you have data points and the corresponding vector of data.

The output  is an interpolation function you can call at any coordinate in a hypercube defined by the bounding box of your input data or the domain explicitely defined by keyword xi

If the user want more control on the grid he needs at least to provide xi, here with option to just provide vectors of coordinates in each direction (only works anyway for this case)  so xi can be a tuple of vectors
    
You can use all keyword parameters of divand

# Input:
* `x`: tuple of arrays of coordinates

* 
* `f`: the value of the function to interpolate at the coordinates x`

# Output:
* `myfunny`: the interpolation function. If you had two dimensional input (i.e. x was a tuple of two coordinates), you can evaluate the interpolation as myfunny(0.1,0.2) for example

"""	
	
	
	
	
	
	
	
	
    maxsize=[1000,100,30,20]
    maxs=10
    if length(x)>4
        @warn "Are you sure ?"
    
        else
        maxs=maxsize[length(x)]
    end
    
    if size(f)[1]!=size(x[1])[1]
        @error "Need same size for coordinates than values"
    end
    # Determine grid resolution 
    nd=size(f)[1]
    
    
    
    NX=min(5*Int(round(10^(log10(nd)/length(x)))),maxs)*ones(Int32,length(x))
    
    # Allocate array
    LX=ones(Float64,length(x))
    
    
    
    
    coords=[]
    if xi==nothing
        
        for i=1:length(x)
            ri=[extrema(x[i])...]
            LX[i]=ri[2]-ri[1]
            # slightly enlarge to allow slight roundings
            ri[2]=ri[2]+0.0001*LX[i]
            ri[1]=ri[1]-0.0001*LX[i]
            LX[i]=ri[2]-ri[1]
            coords=[coords...,range(ri...,NX[i])]
        end
        mask,pmn,xi = DIVAnd_rectdom(coords...)
      else
        
        # No check done just taking the first point...
        for i=1:length(x)
            if isa(xi[i],Vector)
            # Get the coordinates along that direction
             coorxi=deepcopy(xi[i])    
            else
            # Extract grid position in each direction assuming separable coordinates.
             coorxi=deepcopy(xi[i][1:stride(xi[i],i):stride(xi[i],i)*size(xi[i])[i]])
            end
            coords=[coords...,coorxi]
        end
        
        if isa(xi[1],Vector)
            mask,pmn,xi = DIVAnd_rectdom(coords...)
        end
                
        
        if pmn==nothing
            # metric (inverse of the resolution)
            pmn = ndgrid([1 ./ localresolution(coords[i]) for i = 1:length(coords)]...)
        end
        
        
        if mask==nothing
            mask=trues(size(xi[1]))
        end  
    end
    
    if len==nothing
        len=tuple(LX./NX .* 4.0 ...)
    end
    if epsilon2==nothing
        epsilon2=0.1
    end
    
    # Now make the DIVAndrun call
    backg=sum(f)/size(f)[1]
    fi,s=DIVAndrun(mask,pmn,xi,x,f.-backg,len,epsilon2; kwargs...)
    
    # Now initialize interpolations function and return that
     
    return LinearInterpolation(tuple(collect.(coords)...), fi.+backg)
end