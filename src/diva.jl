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
                transform = Anam.notransform(),
                )

    # metadata of grid
    lonr,latr,depthr,TS = xi
    
    # metadata of observations
    lon,lat,depth,time = x

    # correlation length
    lenx,leny,lenz = map(x -> Float64.(x),len)

    # anamorphosis transform
    trans,invtrans = transform
    
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

        # apply the transformation
        value_trans = trans.(value[sel])
        
        # spatial mean of observations
        vm = mean(value_trans)
        va = value_trans - vm
        
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

        # inverse anamorphosis transformation
        fit .= invtrans.(fit)
        
        plotres(timeindex,sel,fit,erri)
        
        divand.writeslice(ncvar, ncvar_relerr, ncvar_Lx,
                          fit, erri, (:,:,:,timeindex))
        
        # write to file
        sync(ds)
    end

    close(ds)
end    



