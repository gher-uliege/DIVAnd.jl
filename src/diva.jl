function diva(dimnames,grid,obscoord,obsvalue,epsilon2,lens,timeaggregation2;
              datadir = joinpath(dirname(@__FILE__),"..","..","divand-example-data"),
              bathname = joinpath(datadir,"Global","Bathymetry","gebco_30sec_16.nc"),
              bathisglobal = true,
              defaultdepth = 0
              )

    # check parameters
    if (length(dimnames) != length(grid)) &&
        (length(dimnames) != length(obscoord))
        
        throw(ArgumentError("the tuples dimnames, grid and obscoord must have the same length"))
    end        
        
    
    for i = 1:length(grid)
        if length(grid[i]) == 0
            throw(ArgumentError("the element $(i) in grid has a length of zero"))
        end

        if length(obscoord[i]) != length(obscoord[1])
            throw(ArgumentError("the element in obscoord must have the same length"))
        end
    end
            
    
    obscoord_aggregated = [obscoord...]
    
    #g = Dict{Symbol,Union{Array{eltype(grid[1]),1},FloatRange{eltype(grid[1])}}}();
    g = Dict{Symbol,Any}();

    for i = 1:length(obscoord)
        if dimnames[i] == "time"
            #@show i
            obscoord_aggregated[i][:] = timeaggregation2(obscoord[i])
            #@show obscoord_aggregated[4][1:10]
            #@show timeaggregation2(obscoord[i][1:4])
        end

        g[Symbol(dimnames[i])] = grid[i]
    end

    #@show obscoord_aggregated[4][1:10]
    mask = 
        if "depth" ∈ Set(dimnames)
            load_mask(bathname,bathisglobal,
                      g[:longitude],
                      g[:latitude],
                      g[:depth])[3]
        else
            # if no depth is provided use defaultdepth 
            load_mask(bathname,bathisglobal,
                      g[:longitude],
                      g[:latitude],
                      defaultdepth)[3]
        end
                
    if "time" ∈ Set(dimnames)
        if  "depth" ∈ Set(dimnames)
            mask = repeat(mask,inner = (1,1,1,length(g[:time])))
        else
            mask = repeat(mask,inner = (1,1,length(g[:time])))
        end
    end

    # size of the analysed array
    sz = ([length(c) for c in grid]...)

    # turns len into arrays (if necessary)
    lens = ([( isa(l,Number) ? fill(l,sz) : l)  for l in lens]...)
    
    fi = SharedArray{eltype(grid[1])}(sz)
    #s = SharedArray(divand.divand_struct,(sz[3],sz[4]))

    #@show size(mask)
    #@show size(grid[1])

    if "longitude" ∈ Set(dimnames) && "latitude" ∈ Set(dimnames)
        xi,yi = ndgrid(g[:longitude],g[:latitude])
        pm,pn = divand_metric(xi,yi)
    end        
    
    
    if dimnames == ("longitude","latitude","depth","time")
        lonr,latr,depthr,timer = grid
        lon,lat,depth,time = obscoord_aggregated
        lenx,leny,lenz,lent = lens


        # begin
        #     @sync @parallel for n = 1:sz[4]
        #         #@show n
        #         for k = 1:sz[3]
        #             sel = (depth .== depthr[k]) & (time .== timer[n])
        #             vm = mean(obsvalue[sel])
        #             va = obsvalue[sel] - vm
                    
        #             #@show sum(sel)
        #             #@show time[1:10]
        #             #@show sum((time .== timer[n]))
        #             #@show vm
        #             #@show va
                                        
        #             fi[:,:,k,n],stmp = divandrun(mask[:,:,k,n],(pm,pn),(xi,yi),(lon[sel],lat[sel]),va,(lenx[:,:,k,n],leny[:,:,k,n]),epsilon2)
        #             fi[:,:,k,n] = fi[:,:,k,n] + vm
        #         end
        #     end
        # end
                     
 
        begin
            @sync @parallel for n = 1:sz[4]
                #@show n
                for k = 1:sz[3]
                    sel = (depth .== depthr[k]) & (time .== timer[n])
                    vm = mean(obsvalue[sel])
                    va = obsvalue[sel] - vm
                    
                    #@show sum(sel)
                    #@show time[1:10]
                    #@show sum((time .== timer[n]))
                    #@show vm
                    #@show va
                                        
                    fi[:,:,k,n],stmp = varanalysis(mask[:,:,:,n],(pm,pn),(xi,yi),(lon[sel],lat[sel]),va,(lenx[:,:,k,n],leny[:,:,k,n]),epsilon2)
                    fi[:,:,k,n] = fi[:,:,k,n] + vm
                end
            end
        end

        #@show fi[:,:,1,1]
        return fi
    elseif dimanmes == ("longitude","latitude","time")
        lenx,leny,lent = lens
        lon,lat,time = obscoord_aggregated
 
        @time begin
            @sync for k = 1:sz[3]
                sel = time .== timer[n]
                vm = mean(obsvalue[sel])
                va = obsvalue[sel] - vm
                    
                #@show sum(sel)
                fi[:,:,k],stmp = divandrun(mask[:,:,k],(pm,pn),(xi,yi),(lon[sel],lat[sel]),va,(lenx[:,:,k],leny[:,:,k]),epsilon2)                
                fi[:,:,k] = fi[:,:,k] + vm
            end
        end
        
        return fi


    end

end


"""Aggregate according to the time `time` (DateTime structure). The function return the aggregated time vector where
all elements have the same value which should be considered togeter. """
aggregation_monthly(time) =  Dates.month.(time)
