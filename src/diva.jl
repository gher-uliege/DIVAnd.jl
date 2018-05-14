"""
    residuals = diva3d(xi,x,value,len,epsilon2,filename,varname)

Create a 3D analysis (or a series of 3D analyses) with DIVAnd using the
observations `value` (vector) at the locations `x` (tuple of vectors) onto
the regular grid defined by the vectors `xi` using the scaled observational error
variance  `epsilon2` and the correlation length `len`. The result will be saved in the
NetCDF file `filename` under the variable `varname`.

## Inputs

*  `xi`: tuple with n elements. Every element represents a coordinate
  of the final grid on which the observations are interpolated

* `x`: tuple with n elements. Every element represents a coordinate of
  the observations

* `value`: value of the observations

* `len`: tuple with n elements. Every element represents the correlation length.
   If `fitcorrlen` is `true`, then `len` can be the empty tuple `()` or a tuple containing
   3 arrays of normalized correlation lengths which will be multiplied by the
   horizontal and vertical correlation lengths.

* `epsilon2`: error variance of the observations (normalized by the error variance of the background field). `epsilon2` can be a scalar (all observations have the same error variance and their errors are decorrelated), a vector (all observations can have a different error variance and their errors are decorrelated) or a matrix (all observations can have a different error variance and their errors can be correlated). If `epsilon2` is a scalar, it is thus the *inverse of the signal-to-noise ratio*.

* `filename`: The output NetCDF filename.

* `varname`: The name of the variable (used in the NetCDF file).

## Optional input arguments:

* `bathname`: path to the NetCDF bathymetry (default ../../divand-example-data/Global/Bathymetry/gebco_30sec_16.nc relative to this source file)
* `bathisglobal`: true (default) is the bathymetry is a global data set
* `plotres`: Call-back routine for plotting ((timeindex,sel,fit,erri) -> nothing)
* `timeorigin`: Time origin (default DateTime(1900,1,1,0,0,0))
* `moddim`: modulo for cyclic dimension (vector with n elements).
     Zero is used for non-cyclic dimensions. Halo points should
     not be included for cyclic dimensions. For example if the first dimension
     is cyclic, then the grid point corresponding to `mask[1,j]` should be
     between `mask[end,1]` (left neighbor) and `mask[2,j]` (right neighbor). The default is [0,0,0],
* `zlevel`: `:surface` (default) for surface analysis and `:floor` for analysis from the bottom floor.
* `ncvarattrib`: dictionary of NetCDF variable attributes.
* `ncglobalattrib`: dictionary of NetCDF global attributes.
* `transform`: Anamorphosis transformation function (default: `Anam.notransform()`).
* `fitcorrlen`: true of the correlation length is determined from the observation (default `false`).
* `fithorz_param`: dictionary with additional optional parameters for `fithorzlen`.
* `fitvert_param`: dictionary with additional optional parameters for `fitvertlen`.
* `distfun`: function to compute the distance (default `(xi,xj) -> divand.distance(xi[2],xi[1],xj[2],xj[1])`).
* `mask`: if different from `nothing`, then this mask overrides land-sea mask based on the bathymetry
(default `nothing`).
* `background`: if different from `nothing`, then this parameter allows one
to load the background from a call-back function
(default `nothing`).
* `background_espilon2_factor`: multiplication for `epsilon2` when computing the background (default 10.).
* `memtofit`: keyword controlling how to cut the domain depending on the memory
    remaining available for inversion. It is not total memory (default 3).
* `niter_e`: Number of iterations to estimate the optimal scale factor of
   `epsilon2` using Desroziers et al. 2005 (doi: 10.1256/qj.05.108). The default
    is 1 (i.e. no optimization is done).

Any additional keywoard arguments understood by `divandgo` can also be used here
(e.g. velocity constrain)

"""

function diva3d(xi,x,value,len,epsilon2,filename,varname;
                datadir = joinpath(dirname(@__FILE__),"..","..","divand-example-data"),
                bathname = joinpath(datadir,"Global","Bathymetry","gebco_30sec_16.nc"),
                bathisglobal = true,
                plotres = (timeindex,sel,fit,erri) -> nothing,
                timeorigin = DateTime(1900,1,1,0,0,0),
                moddim = [0,0,0],
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

    # metadata of grid
    lonr,latr,depthr,TS = xi

    # metadata of observations
    lon,lat,depth,time = x

    # anamorphosis transform
    trans,invtrans = transform

    mask2,(pm,pn,po),(xi,yi,zi) = divand.domain(
        bathname,bathisglobal,lonr,latr,depthr;
        zlevel = zlevel)

    # allow to the user to override the mask
    if mask == nothing
        mask = mask2
    else
        if size(mask) != size(mask2)
            error("expecting a mask of the size $(size(mask)), " *
                  "but got a mask of the size $(size(mask2))")
        end
        # use mask in the following and not mask2
    end

    sz = size(mask)

    # change the depth of the observation
    if zlevel == :floor
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

    # correlation length
    if len == ()
        if !fitcorrlen
            error("if len is () (the empty tuple) then fitcorrlen must be true")
        end

        lenx = ones(sz)
        leny = ones(sz)
        lenz = ones(sz)
    else
        lenx,leny,lenz = map(x -> Float64.(x),len)
    end
    # unscaled copy
    len0 = deepcopy((lenx,leny,lenz))


    # epsilon
    epsilon2 =
        if ndims(epsilon2) == 0
            fill(epsilon2,size(value))
        end

    # days since timeorigin
    timeclim = [Dates.Millisecond(tc - timeorigin).value/(1000*60*60*24) for tc in ctimes(TS)]
    climatologybounds = climatology_bounds(TS)

    # create the NetCDF file
    info("Creating netCDF file")
    ds, ncvar, ncvar_relerr, ncvar_Lx = divand.ncfile(
        filename,(lonr,latr,depthr,timeclim),varname;
        ncvarattrib = ncvarattrib,
        ncglobalattrib = ncglobalattrib,
        climatology_bounds = climatologybounds,
        relerr = true)

    # Prepare background as mean vertical profile and time evolution.
    # Just call divand in two dimensions forgetting x and y ...
    toaverage = [true,true,false]

    # fit fitting
    dbinfo = Dict{Symbol,Any}()

    if fitcorrlen
        kmax = length(depthr)

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

    dbinfo[:factore] = zeros(niter_e,length(TS))

    info("Starting loop on time indices")
    for timeindex = 1:length(TS)

        info("Step $(timeindex) / $(length(TS))")
        info("------------")

        # select observation to be used for the time instance timeindex
        sel = select(TS,timeindex,time)

        if sum(sel) == 0
            warning("no data at $(timeindex)")

            fit = zeros(sz)
            erri = ones(sz)

            divand.writeslice(ncvar, ncvar_relerr, ncvar_Lx,
                              fit, erri, (:,:,:,timeindex))

            # write to file
            sync(ds)

            # next iteration
            continue
        end

        # apply the transformation
        value_trans = trans.(value[sel])


        fbackground,vaa =
            if background == nothing
                info("Computing background profile")
                # spatial mean of observations
                vm = mean(value_trans)
                va = value_trans - vm

                # background profile
                fi,vaa = divand.divand_averaged_bg(
                    mask,(pm,pn,po),(xi,yi,zi),
                    (lon[sel],lat[sel],depth[sel]),va,
                    background_len,
                    epsilon2[sel]*background_epsilon2_factor,
                    toaverage;
                    moddim = moddim)

                fbackground = fi + vm

                fbackground,vaa
            else
                # anamorphosis transform must already be included to the
                # background
                background((lon[sel],lat[sel],depth[sel]),timeindex,value_trans,trans)
            end

        if fitcorrlen
            info("Applying fit of the correlation length")
            # fit correlation length
            lenxy1,infoxy = divand.fithorzlen(
                (lon[sel],lat[sel],depth[sel]),vaa,depthr;
                distfun = distfun,
                fithorz_param...
            )

            lenz1,infoz = divand.fitvertlen(
                (lon[sel],lat[sel],depth[sel]),vaa,depthr;
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
            for k = 1:size(lenx,3)
                for j = 1:size(lenx,2)
                    for i = 1:size(lenx,1)
                        lenx[i,j,k] = len0[1][i,j,k] * lenxy1[k] * pi/180 * EarthRadius
                        leny[i,j,k] = len0[2][i,j,k] * lenxy1[k] * pi/180 * EarthRadius
                        lenz[i,j,k] = len0[2][i,j,k] * lenz1[k]
                    end
                end
            end
        end


        # factore is the total (cumulative) scale factor for
        # espilon2 (Desroziers)

        factore = 1.
        fi2 = zeros(sz)
        erri = zeros(sz)

        for i = 1:niter_e
            info("Estimating the optimal scale factor of epsilon2")
            # error only require at the last iterations
            errortype = (i == niter_e ? :cpme : :none)

            # analysis
            fi2, erri, residuals[sel], qcvalues, scalefactore =
                divand.divandgo(mask,(pm,pn,po),(xi,yi,zi),
                                (lon[sel],lat[sel],depth[sel]),
                                vaa,
                                (lenx,leny,lenz),
                                factore * epsilon2[sel],
                                errortype;
                                moddim = moddim, MEMTOFIT = memtofit, kwargs...)

            factore = scalefactore * factore

            dbinfo[:factore][i,timeindex] = factore
        end

        #fi2,s = divand.varanalysis(mask,(pm,pn,po,pp),(xi,yi,zi,ti),(lon,lat,depth,time2),vaa,(lenx,leny,lenz,lent),epsilon2[sel];                          progress = divand.cgprogress, tol = tol)

        # sum analysis and backgrounds
        fit = fi2 + fbackground

        # inverse anamorphosis transformation
        fit .= invtrans.(fit)

        plotres(timeindex,sel,fit,erri)

        divand.writeslice(ncvar, ncvar_relerr, ncvar_Lx,
                          fit, erri, (:,:,:,timeindex))

        # write to file
        sync(ds)
    end

    close(ds)
    dbinfo[:residuals] = residuals

    return dbinfo
end
