"""
    residuals = diva3d(xi,x,value,len,epsilon2,filename,varname)

Create a 3D analysis (or a series of 3D analyses) with DIVAnd using the
observations `value` (vector) at the locations `x` (tuple of vector) onto
the regular grid defined by the vectors `xi` using the scaled obs. error
variance  `epsilon2` and the correlation length `len`. The result will be saved in the
NetCDF file `filename` under the variable `varname`.

## Inputs

*  `xi`: tuple with n elements. Every element represents a coordinate
  of the final grid on which the observations are interpolated

* `x`: tuple with n elements. Every element represents a coordinate of
  the observations

* `value`: value of the observations

* `len`: tuple with n elements. Every element represents the correlation length.
   If `fitcorrlen`, then `len` can be the emply tuple `()` or a tuple containing
   3 arrays of normalized correlation length which will be multiplied by the
   horizontal and vertical correlation length.

* `epsilon2`: error variance of the observations (normalized by the error variance of the background field). `epsilon2` can be a scalar (all observations have the same error variance and their errors are decorrelated), a vector (all observations can have a different error variance and their errors are decorrelated) or a matrix (all observations can have a different error variance and their errors can be correlated). If `epsilon2` is a scalar, it is thus the *inverse of the signal-to-noise ratio*.

* `filename`: The output NetCDF filename

* `varname`: The name of the variable (used in the NetCDF file)

## Optional input arguments:

* `bathname`: path to the NetCDF bathymetry (default ../../divand-example-data/Global/Bathymetry/gebco_30sec_16.nc relative to this source file)
* `bathisglobal`: true (default) is the bahtymetry is a global data set
* `plotres`: Call-back routine for plotting ((timeindex,sel,fit,erri) -> nothing)
* `timeorigin`: Time origin (default DateTime(1900,1,1,0,0,0))
* `moddim`: modulo for cyclic dimension (vector with n elements).
     Zero is used for non-cyclic dimensions. Halo points should
     not be included for cyclic dimensions. For example if the first dimension
     is cyclic, then the grid point corresponding to `mask[1,j]` should be
     between `mask[end,1]` (left neighbor) and `mask[2,j]` (right neighbor). The default is [0,0,0],
* `zlevel`: :surface (default) for surface analysis and :floor for analysis from the bottom floor
* `ncvarattrib`: Dict of NetCDF variable attributes
* `ncglobalattrib`: Dict of NetCDF global attributes
* `transform`: Anormphosis transformation function (default: `Anam.notransform()`)
* `fitcorrlen`: true of the correlation length is determined from the observation (default `false`)
* `fithorz_param`: is a dictionary with additional optional parameters for `fithorzlen`
* `fitvert_param`: is a dictionary with additional optional parameters for `fitvertlen`
* `distfun`: function to compute the distance (default `(xi,xj) -> divand.distance(xi[2],xi[1],xj[2],xj[1])`)
* `mask`: if different from nothing, then this mask overrride land-sea mask based on the bathymetry (default `nothing`)
* `background`: if different form, then this parameter allows to load the background from a call-back function (default `nothing`)
* `background_espilon2_factor`: multiplication for `epsilon2` when computing the background (default 10.)
* `memtofit`: keyword controlling how to cut the domain depending on the memory
    remaining available for inversion. It is not total memory. (default 3)
* `niter_e`: Number of uterations to estimate the optimal scale factor of 
   `epsilon2` using Desroziers et al. 2005 (doi: 10.1256/qj.05.108). The default
    is 1 (i.e. no optimization is done).

Any additional keywoard arguments understood by divandgo can also be used here 
(e.g. velocity constrain)

"""

function diva3d(xi,x,value,len,epsilon2,filename,varname;
                datadir = joinpath(dirname(@__FILE__),"..","..","divand-example-data"),
                bathname = joinpath(datadir,"Global","Bathymetry","gebco_30sec_16.nc"),
                bathisglobal = true,
                plotres = (timeindex,sel,fit,erri) -> nothing,
                timeorigin = DateTime(1900,1,1,0,0,0),
                moddim = zeros(length(xi)-1),
                zlevel = :surface,
                ncvarattrib = Dict(),
                ncglobalattrib = Dict(),
                transform = Anam.notransform(),
                distfun = (xi,xj) -> divand.distance(xi[2],xi[1],xj[2],xj[1]),
                mask = nothing,
                background = nothing,
                background_epsilon2_factor::Float64 = 10.,
                background_len = (len[1],len[2],4*len[3]),
                fitcorrlen::Bool = false,
                fithorz_param = Dict(),
                fitvert_param = Dict(
                    :distbin => collect([0.:50:400; 500:100:600]),
                    :nmean => 500,
                ),
                memtofit = 3,
                niter_e::Int = 1,
                kwargs...
                )

    # dimension of the analysis
    n = length(xi)
    
    # metadata of grid
    lonr,latr,depthr,TS =
        if length(xi) == 4
            xi
        else
            (xi[1],xi[2],Float64[0.],xi[3])
        end

    # metadata of observations
    lon,lat,depth,time =
        if length(xi) == 4
            x
        else
            (x[1],x[2],Float64[0.],x[3])
        end


    # anamorphosis transform
    trans,invtrans = transform

    # xi[end] is time, which is skipped for the definition of the domain
    mask2,pmn,xyi =
        if length(xi) == 4
            divand.domain(
                bathname,bathisglobal,xi[1],xi[2],xi[3];
                zlevel = zlevel)
        else
            divand.domain(bathname,bathisglobal,xi[1],xi[2])            
        end

    # allow to the user to override the mask
    if mask == nothing
        mask = mask2
    else
        if size(mask) != size(mask2)
            error("expecting a mask of the size $(size(mask2)), " *
                  "but got a mask of the size $(size(mask))")
        end
        # use mask in the following and not mask2
    end

    sz = size(mask)

    # change the depth of the observation
    if (zlevel == :floor) && (n == 4)
        depth = copy(depth)
        info("analysis from the sea floor")

        bxi,byi,bi = divand.load_bath(bathname,bathisglobal,lonr,latr)

        itp = interpolate((bxi,byi), bi, Gridded(Linear()))

        # shift the depth of the observations relative to the ocean floor
        for k = 1:length(depth)
            depth[k] = itp[lon[k],lat[k]] - depth[k]
        end
    end

    # allocate residuals
    residuals = similar(value)
    qcvalues = similar(value)

    # correlation length

    len0 = 
        if len == ()
            if !fitcorrlen
                error("if len is () (the empty tuple) then fitcorrlen must be true")
            end

            # unscaled copy
            ntuple(i -> ones(sz),n-1)
        else
            map(x -> Float64.(x),len)
        end
    # scaling comes later
    len_scaled = deepcopy(len0)
    
    # episilon
    epsilon2 = 
        if ndims(epsilon2) == 0
            fill(epsilon2,size(value))
        end
    
    # days since timeorigin
    timeclim = [Dates.Millisecond(tc - timeorigin).value/(1000*60*60*24) for tc in ctimes(TS)]
    climatologybounds = climatology_bounds(TS)

    # create the NetCDF file
    ds, ncvar, ncvar_relerr, ncvar_Lx =
            divand.ncfile(
                filename,(xi[1:end-1]...,timeclim),varname;
                ncvarattrib = ncvarattrib,
                ncglobalattrib = ncglobalattrib,
                climatology_bounds = climatologybounds,
                relerr = true)    

    # Prepare background as mean vertical profile and time evolution.
    # Just call divand in two dimensions forgetting x and y ...
    toaverage =
        if n == 4
            [true,true,false]
        else
            [true,true]            
        end
    

    # fit fitting
    dbinfo = Dict{Symbol,Any}()

    if fitcorrlen
        kmax = length(depthr)

        if n == 4
            # vertical info
            pmax = length(fitvert_param[:distbin])-1
            dbinfo[:fitvertlen] = Dict{Symbol,Any}(
                :len => zeros(kmax,length(TS)),
                :var0 => zeros(kmax,length(TS)),
                :covar => zeros(pmax,kmax,length(TS)),
                :stdcovar => zeros(pmax,kmax,length(TS)),
                :fitcovar => zeros(pmax,kmax,length(TS)),
                :distx => zeros(pmax)
            )
        end
    end
    
    dbinfo[:factore] = zeros(niter_e,length(TS))
    
    for timeindex = 1:length(TS)
        # select observation to be used for the time instance timeindex
        sel = select(TS,timeindex,time)

        if sum(sel) == 0
            warning("no data at $(timeindex)")

            fit = zeros(sz)
            erri = ones(sz)
            fit[.!mask] = NaN
            erri[.!mask] = NaN

            if n == 4
                divand.writeslice(ncvar, ncvar_relerr, ncvar_Lx,
                                  fit, erri, (:,:,:,timeindex))
            else
                divand.writeslice(ncvar, ncvar_relerr, ncvar_Lx,
                                  fit, erri, (:,:,timeindex))                
            end

            # write to file
            sync(ds)

            # next iteration
            continue
        end

        # apply the transformation
        value_trans = trans.(value[sel])


        # obs. coordinate matching selection
        xsel = map(xc -> xc[sel],x[1:end-1])

        fbackground,vaa =
            if background == nothing
                # spatial mean of observations
                vm = mean(value_trans)
                va = value_trans - vm

                # background profile
                fi,vaa = divand.divand_averaged_bg(
                    mask,pmn,xyi,
                    xsel,
                    va,
                    background_len,
                    epsilon2[sel]*background_epsilon2_factor,
                    toaverage;
                    moddim = moddim)

                fbackground = fi + vm

                fbackground,vaa
            else
                # anamorphosis transform must already be included to the
                # background
                background(xsel,timeindex,value_trans,trans)
            end

        if fitcorrlen
            # fit correlation length
            lenxy1,infoxy = divand.fithorzlen(
                xsel,vaa,depthr;
                distfun = distfun,
                fithorz_param...
            )            

            if n == 3
                # propagate
                for j = 1:sz[2]
                    for i = 1:sz[1]
                        len_scaled[1][i,j] = len0[1][i,j] * lenxy1[1] * pi/180 * EarthRadius
                        len_scaled[2][i,j] = len0[2][i,j] * lenxy1[1] * pi/180 * EarthRadius
                    end
                end
            end

            if n == 4
                lenz1,infoz = divand.fitvertlen(
                    xsel,vaa,depthr;
                    distfun = distfun,
                    fitvert_param...
                )

                dbinfo[:fitvertlen][:len][:,timeindex] = infoz[:len]
                dbinfo[:fitvertlen][:var0][:,timeindex] = infoz[:var0]
                dbinfo[:fitvertlen][:distx][:] = infoz[:distx]
                for key in [:covar,:fitcovar,:stdcovar]
                    dbinfo[:fitvertlen][key][:,:,timeindex] = infoz[key]
                end

                
                # propagate
                for k = 1:sz[3]
                    for j = 1:sz[2]
                        for i = 1:sz[1]
                            len_scaled[1][i,j,k] = len0[1][i,j,k] * lenxy1[k] * pi/180 * EarthRadius
                            len_scaled[2][i,j,k] = len0[2][i,j,k] * lenxy1[k] * pi/180 * EarthRadius
                            len_scaled[3][i,j,k] = len0[3][i,j,k] * lenz1[k]
                        end
                    end
                end
            end
        end


        # factore is the total (cumulative) scale factor for
        # espilon2 (Desroziers)

        factore = 1.
        fi2 = zeros(sz)
        erri = zeros(sz)

        kwargs_without_qcm = [(p,v) for (p,v) in kwargs if p !== :QCMETHOD]
        
        for i = 1:niter_e
            # error and QCMETHOD is only required at the last iterations
            errortype,kwargs2 =
                if i == niter_e
                    :cpme,kwargs
                else
                    :none,kwargs_without_qcm
                end

            # analysis
            fi2, erri, residuals[sel], qcdata, scalefactore =
                divand.divandgo(mask,pmn,xyi,
                                xsel,
                                vaa,
                                len_scaled,
                                factore * epsilon2[sel],
                                errortype;
                                moddim = moddim, MEMTOFIT = memtofit, kwargs2...)

            factore = scalefactore * factore

            dbinfo[:factore][i,timeindex] = factore

            if qcdata != ()
                qcvalues[sel] = qcdata
            end
        end

        #fi2,s = divand.varanalysis(mask,(pm,pn,po,pp),(xi,yi,zi,ti),(lon,lat,depth,time2),vaa,(lenx,leny,lenz,lent),epsilon2[sel];                          progress = divand.cgprogress, tol = tol)

        # sum analysis and backgrounds
        fit = fi2 + fbackground

        # inverse anamorphosis transformation
        fit .= invtrans.(fit)

        plotres(timeindex,sel,fit,erri)

        fit[.!mask] = NaN
        erri[.!mask] = NaN
        if n == 4
            divand.writeslice(ncvar, ncvar_relerr, ncvar_Lx,
                              fit, erri, (:,:,:,timeindex))
        else
            divand.writeslice(ncvar, ncvar_relerr, ncvar_Lx,
                              fit, erri, (:,:,timeindex))
        end

        # write to file
        sync(ds)
    end

    close(ds)
    dbinfo[:residuals] = residuals
    if haskey(Dict(kwargs), :QCMETHOD)
       dbinfo[:qcvalues] = qcvalues
    end

    return dbinfo
end
