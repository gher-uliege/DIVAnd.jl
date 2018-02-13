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


"""
zlevel :surface or :floor
"""
function diva3d(xi,x,value,epsilon2,len,filename,varname;
                datadir = joinpath(dirname(@__FILE__),"..","..","divand-example-data"),
                bathname = joinpath(datadir,"Global","Bathymetry","gebco_30sec_16.nc"),
                bathisglobal = true,
                plotres = (timeindex,sel,fit,erri) -> nothing,
                timeorigin = DateTime(1900,1,1,0,0,0),
                moddim = [0,0,0],
                zlevel = :surface,
                ncvarattrib = Dict(), ncglobalattrib = Dict(),
                )

    # metadata of grid
    lonr,latr,depthr,TS = xi
    
    # metadata of observations
    lon,lat,depth,time = x

    # correlation length
    lenx,leny,lenz = map(x -> Float64.(x),len)
    
    mask,(pm,pn,po),(xi,yi,zi) = divand.domain(
        bathname,bathisglobal,lonr,latr,depthr;
        zlevel = zlevel)

    sz = size(mask)

    if zlevel == :floor       
        depth = copy(depth)
        @show "analysis from the sea floor"
        
        bxi,byi,bi = divand.load_bath(bathname,bathisglobal,lonr,latr)

        itp = interpolate((bxi,byi), bi, Gridded(Linear()))
        
        # shift the depth of the observations relative to the ocean floor
        for k = 1:length(depth)
            depth[k] = -itp[lon[k],lat[k]] - depth[k]
        end
    end

    
    # days since timeorigin
    timeclim = [Dates.Millisecond(tc - timeorigin).value/(1000*60*60*24) for tc in ctimes(TS)]
    
    # create the NetCDF file
    ds, ncvar, ncvar_relerr, ncvar_Lx = divand.ncfile(
        filename,(lonr,latr,depthr,timeclim),varname;
        ncvarattrib = ncvarattrib,
        ncglobalattrib = ncglobalattrib,
        relerr = true)

    # Prepare background as mean vertical profile and time evolution.
    # Just call divand in two dimensions forgetting x and y ...
    toaverage = [true,true,false]   
    
    for timeindex = 1:length(TS)
        # select observation to be used for the time instance timeindex
        sel = select(TS,timeindex,time)

        # spatial mean of observations
        vm = mean(value[sel])
        va = value[sel]-vm
        
        # background profile
        fi,vaa = divand.divand_averaged_bg(mask,(pm,pn,po),(xi,yi,zi),
                                           (lon[sel],lat[sel],depth[sel]),va,
                                           (lenx,leny,4*lenz),epsilon2*10,toaverage;
                                           moddim = moddim)
        
        # analysis
        fi2,erri = divand.divandgo(mask,(pm,pn,po),(xi,yi,zi),
                                   (lon[sel],lat[sel],depth[sel]),vaa,
                                   (lenx,leny,lenz),epsilon2,:cpme;
                                   moddim = moddim, MEMTOFIT = 3.0)
    
        #fi2,s = divand.varanalysis(mask,(pm,pn,po,pp),(xi,yi,zi,ti),(lon,lat,depth,time2),vaa,(lenx,leny,lenz,lent),epsilon2;                          progress = divand.cgprogress, tol = tol)   
        
        # sum analysis and backgrounds
        fit = fi2 + fi + vm
        
        plotres(timeindex,sel,fit,erri)
        
        divand.writeslice(ncvar, ncvar_relerr, ncvar_Lx,
                          fit, erri, (:,:,:,timeindex))
        
        # write to file
        sync(ds)
    end

    close(ds)
end    



