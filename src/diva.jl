"""
    dbinfo = diva3d(xi,x,value,len,epsilon2,filename,varname)

Create a 3D analysis (or a series of 3D analysis) with DIVAnd using the
observations `value` (vector) at the locations `x` (tuple of vectors) onto
the regular grid defined by the vectors `xi` using the scaled observational error
variance `epsilon2` and the correlation length `len`. The result will be saved in the
netCDF file `filename` under the variable `varname`.

## Inputs

*  `xi`: tuple with n elements. Every element represents a coordinate
  of the final grid on which the observations are interpolated

* `x`: tuple with n elements. Every element represents a coordinate of
  the observations

* `value`: value of the observations

* `len`: tuple with n elements. Every element represents the correlation length.
   If `fitcorrlen` is `false` (default), the correlation length should be expressed in meters.
   If `fitcorrlen` is `true`, then `len` can be the empty tuple `()` or a tuple containing
   3 arrays of normalized correlation lengths which will be multiplied by the
   horizontal and vertical correlation lengths.

* `epsilon2`: error variance of the observations (normalized by the error variance of the background field).
`epsilon2` can be a scalar (all observations have the same error variance and their errors are decorrelated), a vector (all observations can have a different error variance and their errors are decorrelated) or a matrix (all observations can have a different error variance and their errors can be correlated). If `epsilon2` is a scalar, it is thus the *inverse of the signal-to-noise ratio*.

* `filename`: The output netCDF filename.

* `varname`: The name of the variable (used in the netCDF file).

## Optional input arguments:

* `bathname`: path to the netCDF bathymetry (default ../../DIVAnd-example-data/Global/Bathymetry/gebco_30sec_16.nc relative to this source file)
* `bathisglobal`: true (default) is the bathymetry is a global data set
* `plotres`: Call-back routine for plotting ((timeindex,sel,fit,erri) -> nothing)
* `timeorigin`: Time origin (default DateTime(1900,1,1,0,0,0))
* `moddim`: modulo for cyclic dimension (vector with n elements).
     Zero is used for non-cyclic dimensions. Halo points should
     not be included for cyclic dimensions. For example if the first dimension
     is cyclic, then the grid point corresponding to `mask[1,j]` should be
     between `mask[end,1]` (left neighbor) and `mask[2,j]` (right neighbor). The default is [0,0,0],
* `zlevel`: `:surface` (default) for surface analysis and `:floor` for analysis from the bottom floor.
* `ncvarattrib`: dictionary of netCDF variable attributes.
* `ncglobalattrib`: dictionary of netCDF global attributes.
* `transform`: Anamorphosis transformation function (default: `Anam.notransform()`).
* `fitcorrlen`: true if the correlation length is determined from the observation (default `false`).
     Note that the parameter `len` is interpreted differently when `fitcorrlen` is set to `true`.
* `fithorzcorrlen`: true if the horizontal correlation length is determined from the observation (default: the value of `fitcorrlen`)
     Note that the parameter `len` is interpreted differently when `fithorzcorrlen` is set to `true`.
* `fitvertcorrlen`: true if the vertical correlation length is determined from the observation (default: the value of `fitcorrlen`)
     Note that the parameter `len` is interpreted differently when `fitvertcorrlen` is set to `true`.
* `fithorz_param`: dictionary with additional optional parameters for `fithorzlen`, for example: `Dict(:smoothz => 200., :searchz => 50.)`.
* `fitvert_param`: dictionary with additional optional parameters for `fitvertlen`.
* `distfun`: function to compute the distance (default `(xi,xj) -> DIVAnd.distance(xi[2],xi[1],xj[2],xj[1])`).
* `mask`: if different from `nothing`, then this mask overrides land-sea mask based on the bathymetry
(default `nothing`).
* `background`: if different from `nothing`, then this parameter allows one
to load the background from a call-back function (default `nothing`). The call-back functions has the parameters
`(x,n,trans_value,trans)` where `x` represent the position of the observations, `n` the time index, `trans_value`, the observations
(possibly transformed) and `trans` the transformation function. The output of this function is the
gridded background field and the observations minus the background field.
* `background_espilon2_factor`: multiplication for `epsilon2` when computing a
   vertical profile as a background estimate (default 10.). This parameter is not used
   when the parameter `background` or `background_lenz` is provided.
* `background_lenz`: vertical correlation for background computation (default 20 m). This parameter is not used
   when the parameter `background` is provided.
* `background_len`: deprecated option replaced by `background_lenz`.
* `memtofit`: keyword controlling how to cut the domain depending on the memory
    remaining available for inversion. It is not the total memory (default 3). Use a large value (e.g. 100) to force the
    usage for the more efficient direct solver if you are not limited by the amount of RAM memory.
* `minfield`: if the analysed field is below `minfield`, its value is replace by `minfield` (default -Inf, i.e. no substitution is done).
* `maxfield`: if the analysed field is above `maxfield`, its value is replace by `maxfield` (default +Inf, i.e. no substitution is done).
* `saveindex`: controls if just a subset of the analysis should be saved to
    the netCDF file. Per default, `saveindex` is `(:,:,:)` (corresponding to
    longitude, latitude and depth indices) meaning that everything is saved.
    If however, for example the first layer should not be saved then `saveindex`
    should be `(:,:,2:length(depthr))` where `depthr` is the 3rd element of `xi`.
* `niter_e`: Number of iterations to estimate the optimal scale factor of
   `epsilon2` using Desroziers et al. 2005 (doi: 10.1256/qj.05.108). The default
    is 1 (i.e. no optimization is done).
* `coeff_derivative2` (vector of 3 floats): for every dimension where this value is non-zero, an additional term is added to the cost function penalizing the second derivative. A typical value of this parameter is `[0.,0.,1e-8]`.

Any additional keywoard arguments understood by `DIVAndgo` can also be used here
(e.g. velocity constrain)


The output is a dictionary with the followings keys:

* `:residuals`: the difference between the analysis (interpolated linearly to the
location of the observation) and the observations. The
residual is NaN if the observations are not within the domain as defined by
the mask and the coordinates of the observations `x`.
* `:qcvalues`: quality control scores (if activated)

!!! note

    At all vertical levels, there should at least one sea point.
"""
function diva3d(xi,x,value,len,epsilon2,filename,varname;
                datadir = joinpath(dirname(@__FILE__),"..","..","DIVAnd-example-data"),
                bathname = joinpath(datadir,"Global","Bathymetry","gebco_30sec_16.nc"),
                bathisglobal = true,
                plotres = (timeindex,sel,fit,erri) -> nothing,
                timeorigin = DateTime(1900,1,1,0,0,0),
                moddim = zeros(length(xi)-1),
                zlevel = :surface,
                ncvarattrib = Dict(),
                ncglobalattrib = Dict(),
                transform = Anam.notransform(),
                distfun = distfun_m,
                mask = nothing,
                background = nothing,
                background_epsilon2_factor::Float64 = 10.,
                background_lenz = nothing, # m
                background_len = nothing,
                background_lenz_factor = 4,
                fitcorrlen::Bool = false,
                fithorzcorrlen::Bool = fitcorrlen,
                fitvertcorrlen::Bool = fitcorrlen,
                fithorz_param = Dict(),
                fitvert_param = Dict(),
                memtofit = 3,
                niter_e::Int = 1,
                minfield::Number = -Inf,
                maxfield::Number = Inf,
                surfextend = false,
                kwargs...
                )

    # dimension of the analysis
    n = length(xi)

    # save everything per default
    saveindex = ntuple(i -> :, length(xi)-1)
    background_ext = background
    plotres_ext = plotres

    if surfextend && length(len) == 3
        if all(len[3] .== 0)
            error("surfextend should not be used when the vertical correlation length is set to zero")
        end
    end

    # vertical extension at surface
    if surfextend
        if length(xi) == 3
            error("surfextend can only be true for 3d analyses")
        end
        if mask != nothing
            mask = cat(mask[:,:,1],mask,dims = 3)
        end

        dz = xi[3][2]-xi[3][1]
        xi = ntuple(i -> (i == 3 ? vcat(xi[3][1]-dz,xi[3]) : xi[i] ), length(xi))

        len = ntuple(i ->
                     if typeof(len[i]) <: Array
                     cat(len[i][:,:,1],len[i],dims = 3)
                     else
                     len[i]
                     end, length(len))
        saveindex = (:,:,2:length(xi[3]))

        if background_ext != nothing
            background_ext =
                (args...; kwargs...) -> begin
                    fbackground,vaa = background(args...; kwargs...)
                    fbackground = cat(fbackground[:,:,1],fbackground,dims = 3)
                    return fbackground,vaa
                end
        end

        #TODO extend background_len if present

        # skip first layer when plotting
        plotres_ext(timeindex,sel,fit,erri) = plotres(timeindex,sel,fit[:,:,2:end],erri[:,:,2:end])
    end

    # metadata of grid
    lonr,latr,depthr,TS =
        if length(xi) == 4
            xi
        else
            (xi[1],xi[2],Float64[0.],xi[3])
        end

    checkdepth(depthr)

    if niter_e < 1
        error("niter_e (currently $niter_e) should be larger than 1")
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
            DIVAnd.domain(
                bathname,bathisglobal,xi[1],xi[2],xi[3];
                zlevel = zlevel)
        else
            DIVAnd.domain(bathname,bathisglobal,xi[1],xi[2])
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

    @static if VERSION >= v"0.7.0-beta.0"
        if  any(sum(mask,dims = [1,2]) .== 0)
            minval,minindex = findmin(sum(mask,dims = [1,2])[:])
            error("some slices completely masked: k = $(minindex)")
        end
    end

    # vertical extension at surface
    if surfextend
        if length(xi) == 3
            error("surfextend can only be true for 3d analyses")
        end
        mask[:,:,1] = mask[:,:,2]
    end

    sz = size(mask)

    # change the depth of the observation
    if (zlevel == :floor) && (n == 4)
        depth = copy(depth)
        @info "analysis from the sea floor"

        bxi,byi,bi = DIVAnd.load_bath(bathname,bathisglobal,lonr,latr)

        itp = interpolate((bxi,byi), bi, Gridded(Linear()))

        # shift the depth of the observations relative to the ocean floor
        for k = 1:length(depth)
            depth[k] = itp[lon[k],lat[k]] - depth[k]
        end
    end

    # allocate residuals and qcvalues
    residuals = similar(value)
    qcvalues = similar(value)
    residuals .= NaN
    qcvalues .= NaN

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

    if background_len == nothing
        if (background_lenz == nothing) && (len == ())
            # default value if nothing is provided
            background_lenz = 20 # m
        end

        # 4D analysis (lon,lat,depth,time)
        if n == 4
            if len !== ()
                if background_lenz !== nothing
                    background_len = (len[1], len[2], background_lenz)
                else
                    background_len = (len[1], len[2], background_lenz_factor*len[3])
                end
            end
        # 3D analysis (lon,lat,time)
        else
            if background_lenz !== nothing
                background_len = (1., # unused
                                  1.) # unused
            end
            if len !== ()
                background_len = (len[1],len[2])
            end
        end
    end


    # epsilon: change to vector if scalar provided (same value for all the points)
    if ndims(epsilon2) == 0
        epsilon2 = fill(epsilon2,size(value))
    end

    # fit fitting
    dbinfo = Dict{Symbol,Any}()

    # days since timeorigin
    timeclim = [Dates.Millisecond(tc - timeorigin).value/(1000*60*60*24) for tc in ctimes(TS)]
    climatologybounds = climatology_bounds(TS)

    # create the NetCDF file
    # make sure that the file is closed even if there is an error
    @info "Creating netCDF file"
    Dataset(filename,"c") do ds
        ncvar, ncvar_relerr, ncvar_Lx =
            DIVAnd.ncfile(
                ds,filename,(xi[1:end-1]...,ctimes(TS)),varname;
                ncvarattrib = ncvarattrib,
                ncglobalattrib = ncglobalattrib,
                climatology_bounds = climatologybounds,
                saveindex = saveindex,
                relerr = true)

        # Prepare background as mean vertical profile and time evolution.
        # Just call DIVAnd in two dimensions forgetting x and y ...
        toaverage =
            if n == 4
                [true,true,false]
            else
                [true,true]
            end

        if fitvertcorrlen
            kmax = length(depthr)

            if n == 4
                # vertical info
                dbinfo[:fitvertlen] = Dict{Symbol,Any}(
                    :len => zeros(kmax,length(TS)),
                    :lenf => zeros(kmax,length(TS)),
                    :var0 => zeros(kmax,length(TS)),
                    :fitinfos => Array{Dict{Symbol,Any},2}(undef,kmax,length(TS))
                )
            end
        end

        dbinfo[:factore] = zeros(niter_e,length(TS))

        for timeindex = 1:length(TS)
            @info "Time step $(timeindex) / $(length(TS))"
            # select observation to be used for the time instance timeindex
            sel = select(TS,timeindex,time)

            if sum(sel) == 0
                @warn "no data at $(timeindex)"

                fit = zeros(sz)
                erri = ones(sz)
                fit[.!mask] .= NaN
                erri[.!mask] .= NaN

                if n == 4
                    DIVAnd.writeslice(ncvar, ncvar_relerr, ncvar_Lx,
                                      fit, erri, (:,:,:,timeindex),
                                      saveindex = saveindex,
                                      )
                else
                    DIVAnd.writeslice(ncvar, ncvar_relerr, ncvar_Lx,
                                      fit, erri, (:,:,timeindex),
                                      saveindex = saveindex,
                                      )
                end

                # write to file
                sync(ds)

                # next iteration
                continue
            end

            # apply the transformation
            value_trans = trans.(value[sel])
            # some transformation (e.g. log) can produce -Inf, set these to NaN
            value_trans[.!isfinite.(value_trans)] .= NaN

            # obs. coordinate matching selection
            xsel = map(xc -> xc[sel],x[1:end-1])

            # @info "Computing background profile"
            fbackground,vaa =
                if background == nothing
                    # spatial mean of observations
                    vm = mean(value_trans[isfinite.(value_trans)])
                    va = value_trans .- vm

                    #@show background_len[3][1,1,:], vm
                    #JLD2.@save "/tmp/test.jld2" background_len mask pmn xyi xsel va epsilon2 sel background_epsilon2_factor toaverage  moddim vm
                    #@show "saving"
                    # background profile

                    if n == 4
                       if all(background_len[3] .== 0)
                           @debug "DIVAnd_averaged_bg use twice the resolution as vertical correlation"
                           background_len[3] .= 2 ./ pmn[3];
                       end
                    end

                    fi,vaa = DIVAnd.DIVAnd_averaged_bg(
                        mask,pmn,xyi,
                        xsel,
                        va,
                        background_len,
                        epsilon2[sel]*background_epsilon2_factor,
                        toaverage;
                        moddim = moddim)

                    fbackground = fi .+ vm
                    @debug "fbackground: $(fbackground[1,1,:])"
                    #@show size(fbackground),fbackground[1,1,end]
                    #dbinfo[:background] = fbackground
                    fbackground,vaa
                else
                    # anamorphosis transform must already be included to the
                    # background
                    background_ext(xsel,timeindex,value_trans,trans;
                               selection = sel, obstime = time)
                end

            # the background is not computed for points outside of the domain
            # exclude those in vaa and xsel
            selbackground = isfinite.(vaa)
            vaa = vaa[selbackground]
            xsel = map(xc -> xc[selbackground],xsel)
            # unselect the data points out of the domain
            view(sel,sel)[.!selbackground] .= false

            if fithorzcorrlen
                # @info "Applying fit of the correlation length"
                # fit correlation length
                fithorz_param_sel = Dict{Symbol,Any}(fithorz_param)
                fithorz_param_sel[:epsilon2] = get(fithorz_param,:epsilon2,epsilon2)[sel]

                lenxy1,infoxy = DIVAnd.fithorzlen(
                    xsel,vaa,depthr;
                    distfun = distfun,
                    fithorz_param_sel...,
                )

                if n == 3
                    # propagate
                    for j = 1:sz[2]
                        for i = 1:sz[1]
                            len_scaled[1][i,j] = len0[1][i,j] * lenxy1[1]
                            len_scaled[2][i,j] = len0[2][i,j] * lenxy1[1]
                        end
                    end
                else
                    for k = 1:sz[3]
                        for j = 1:sz[2]
                            for i = 1:sz[1]
                                len_scaled[1][i,j,k] = len0[1][i,j,k] * lenxy1[k]
                                len_scaled[2][i,j,k] = len0[2][i,j,k] * lenxy1[k]
                            end
                        end
                    end
                end
            end

            if (fitvertcorrlen) && (n == 4)
                fitvert_param_sel = Dict{Symbol,Any}(fitvert_param)
                fitvert_param_sel[:epsilon2] = get(fitvert_param,:epsilon2,epsilon2)[sel]

                lenz1,infoz = DIVAnd.fitvertlen(
                    xsel,vaa,depthr;
                    distfun = distfun,
                    fitvert_param_sel...
                )

                dbinfo[:fitvertlen][:lenf][:,timeindex] = lenz1
                dbinfo[:fitvertlen][:len][:,timeindex] = infoz[:len]
                dbinfo[:fitvertlen][:var0][:,timeindex] = infoz[:var0]
                dbinfo[:fitvertlen][:fitinfos][:,timeindex] = infoz[:fitinfos]

                # propagate
                for k = 1:sz[3]
                    for j = 1:sz[2]
                        for i = 1:sz[1]
                            len_scaled[3][i,j,k] = len0[3][i,j,k] * lenz1[k]
                        end
                    end
                end
            end

            for i = 1:n-1
                @info "scaled correlation length (min,max) in dimension $i: $(extrema(len_scaled[i]))"
            end

            # factore is the total (cumulative) scale factor for
            # espilon2 (Desroziers)

            factore = 1.
            fi2 = zeros(sz)
            erri = zeros(sz)

            kwargs_without_qcm = [(p,v) for (p,v) in kwargs if p !== :QCMETHOD]

            for i = 1:niter_e
                #@info "Estimating the optimal scale factor of epsilon2"
                # error and QCMETHOD is only required at the last iterations
                errortype,kwargs2 =
                    if i == niter_e
                        :cpme,kwargs
                    else
                        :none,kwargs_without_qcm
                    end

                # check the resolution
                checkresolution(mask,pmn,len_scaled)

                # analysis
                fi2, erri, residual, qcdata, scalefactore =
                    DIVAnd.DIVAndgo(mask,pmn,xyi,
                                    xsel,
                                    vaa,
                                    len_scaled,
                                    factore * epsilon2[sel],
                                    errortype;
                                    moddim = moddim, MEMTOFIT = memtofit, kwargs2...)

                residuals[sel] = residual
                factore = scalefactore * factore

                dbinfo[:factore][i,timeindex] = factore

                if qcdata != ()
                    qcvalues[sel] = qcdata
                end
            end

            # sum analysis and backgrounds
            fit = fi2 + fbackground
            #@show niter_e
            #JLD2.@save "/tmp/test_fi2.jld" fi2 fit fbackground

            # inverse anamorphosis transformation
            fit .= invtrans.(fit)

            # apply range check
            if minfield != -Inf
                fit[fit .< minfield] .= minfield
            end
            if maxfield != Inf
                fit[fit .> maxfield] .= maxfield
            end

            #JLD2.@save "/tmp/test_fit.jld2" fit fbackground

            plotres_ext(timeindex,isfinite.(residuals) .& sel,fit,erri)

            fit[.!mask] .= NaN
            erri[.!mask] .= NaN

            if n == 4
                DIVAnd.writeslice(ncvar, ncvar_relerr, ncvar_Lx,
                                  fit, erri, (:,:,:,timeindex),
                                  saveindex = saveindex,
                                  )
            else
                DIVAnd.writeslice(ncvar, ncvar_relerr, ncvar_Lx,
                                  fit, erri, (:,:,timeindex),
                                  saveindex = saveindex,
                                  )
            end

            # write to file
            sync(ds)

        end # time loop

    end # implictly closing the NetCDF file

    dbinfo[:residuals] = residuals
    dbinfo[:used] = .!isnan.(residuals)
    if haskey(Dict(kwargs), :QCMETHOD)
       dbinfo[:qcvalues] = qcvalues
    end
    dbinfo[:mask] = mask

    return dbinfo
end
