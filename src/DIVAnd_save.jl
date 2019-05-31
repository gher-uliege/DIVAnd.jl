function DIVAnd_save(filename,mask::AbstractArray{Bool,N},varname,fi) where N

    sz = size(mask)

    Dataset(filename,"c") do ds
        # Dimensions
        ds.dim["time"] = sz[4]
        ds.dim["depth"] = sz[3]
        ds.dim["lat"] = sz[2]
        ds.dim["lon"] = sz[1]

        ncvar = defVar(ds, varname, Float32, ("lon", "lat", "depth", "time"),
                       checksum = checksum)
        ncvar[:,:,:,:] = fi
    end

    return nothing
end


"""
    DIVAnd_save(ds,filename,xyi,fi,varname;
                      ncvarattrib = Dict(), ncglobalattrib = Dict(), ...)

Save the result of the analysis in a NetCDF file .

# Input arguments

* `ds`: the NetCDF dataset
* `filename`: the name of the NetCDF file
* `mask`: binary mask delimiting the domain. true is inside and false outside. For oceanographic application, this is the land-sea mask where sea is true and land is false.
*  `xyi`: tuple with n elements. Every element represents a coordinate
  of the final grid on which the observations are interpolated
*  `fi`: the analysed field
* `varname`: the name of the NetCDF variable

# Optional arguments:
  * `ncglobalattrib`: a dictionary with the global attributes
  * `ncvarattrib`: a dictionary with the variable attributes
  * `relerr`: relative error
  * `timeorigin`: time origin for the time units attribute (default is 1900-01-01 00:00:00)
"""
function ncfile(ds,filename,xyi,varname;
                ncvarattrib = Dict(), ncglobalattrib = Dict(),
                thresholds = [("L1",0.3),("L2",0.5)],
                deflatelevel = 5,
                chunksizes = [100,100,1,1][1:length(xyi)],
                type_save = Float32,
                timeorigin = DateTime(1900,1,1,0,0,0),
                checksum = :fletcher32,
                saveindex = ntuple(i -> :,length(xyi)-1),
                kwargs...)

    function defnD(ds,varname,dims,ncvarattrib)
        local ncvar = defVar(ds,varname, type_save, dims;
                             deflatelevel = deflatelevel,
                             checksum = checksum,
                             chunksizes = chunksizes[1:length(dims)])

        for (k,v) in ncvarattrib
            ncvar.attrib[k] = v
        end

        ncvar.attrib["_FillValue"] = type_save(fillval)
        ncvar.attrib["missing_value"] = type_save(fillval)

        return ncvar
    end

    def4D(ds,varname,ncvarattrib) = defnD(ds,varname,dimnames,ncvarattrib)

    # unused
    # def3D(ds,varname,ncvarattrib) = defnD(ds,varname,("lon", "lat", "time"),ncvarattrib)

    kw = Dict(kwargs)

    # size
    sz = length.(xyi)

    # chunksizes should not exceed the size of fi
    chunksizes = min.(chunksizes,collect(sz))

    # index of depth and time (-1 mean no depth/time dimension)
    idepth = -1
    itime = -1
    dimnames = ("lon","lat")

    if length(xyi) == 4
        idepth = 3
        itime = 4
        dimnames = ("lon","lat","depth","time")
    elseif length(xyi) == 3
        if eltype(xyi[3]) <: DateTime
            itime = 3
            dimnames = ("lon","lat","time")
        else
            idepth = 3
            dimnames = ("lon","lat","depth")
        end
    end

    units = get(ncvarattrib,"units","")
    longname = get(ncvarattrib,"long_name","")
    validmin = get(ncvarattrib,"valid_min","")
    validmax = get(ncvarattrib,"valid_max",0.)

    fillval = NC_FILL_FLOAT

    # Dimensions

    ds.dim["lon"] = length((1:sz[1])[saveindex[1]])
    ds.dim["lat"] = length((1:sz[2])[saveindex[2]])

    if idepth != -1
        ds.dim["depth"] = length((1:sz[idepth])[saveindex[3]])
    end

    if itime != -1
        ds.dim["time"] = sz[itime]
    end

    ds.dim["nv"] = 2

    # Declare variables

    nclon = defVar(ds,"lon", Float64, ("lon",), checksum = checksum)
    nclon.attrib["units"] = "degrees_east"
    nclon.attrib["standard_name"] = "longitude"
    nclon.attrib["long_name"] = "longitude"


    nclat = defVar(ds,"lat", Float64, ("lat",), checksum = checksum)
    nclat.attrib["units"] = "degrees_north"
    nclat.attrib["standard_name"] = "latitude"
    nclat.attrib["long_name"] = "latitude"


    if idepth != -1
        ncdepth = defVar(ds,"depth", Float64, ("depth",), checksum = checksum)
        ncdepth.attrib["units"] = "meters"
        ncdepth.attrib["positive"] = "down"
        ncdepth.attrib["standard_name"] = "depth"
        ncdepth.attrib["long_name"] = "depth below sea level"
    end

    if itime != -1
        nctime = defVar(ds,"time", Float64, ("time",), checksum = checksum)
        nctime.attrib["units"] = "days since " *
            Dates.format(timeorigin,"yyyy-mm-dd HH:MM:SS")
        nctime.attrib["standard_name"] = "time"
        nctime.attrib["long_name"] = "time"
        nctime.attrib["calendar"] = "standard"

        if haskey(kw,:climatology_bounds)
            nctime.attrib["climatology"] = "climatology_bounds"
            ncclimatology_bounds = defVar(ds,"climatology_bounds", Float64,
                                          ("nv", "time"), checksum = checksum)
            ncclimatology_bounds.attrib["units"] = "days since " *
                Dates.format(timeorigin,"yyyy-mm-dd HH:MM:SS")
        end
    end


    # ncCLfield = defVar(ds,"CLfield", type_save, ("lon", "lat", "depth", "time"))
    # ncCLfield.attrib["long_name"] = "Correlation length field"
    # ncCLfield.attrib["valid_min"] = type_save(0.0)
    # ncCLfield.attrib["valid_max"] = type_save(0.76)
    # ncCLfield.attrib["_FillValue"] = type_save(fillval)
    # ncCLfield.attrib["missing_value"] = type_save(fillval)

    # ncCORRLEN = defVar(ds,"CORRLEN", type_save, ("depth", "time"))
    # ncCORRLEN.attrib["long_name"] = "Correlation Length"
    # ncCORRLEN.attrib["units"] = "degrees_north"

    # ncSNR = defVar(ds,"SNR", type_save, ("depth", "time"))
    # ncSNR.attrib["long_name"] = "Signal to Noise"

    # ncVARBACK = defVar(ds,"VARBACK", type_save, ("depth", "time"))
    # ncVARBACK.attrib["long_name"] = "Background Field Variance"
    # ncVARBACK.attrib["units"] = "umol/l^2"


    ncvar = def4D(ds,varname,ncvarattrib)
    ncvar.attrib["long_name"] = "$(longname)"
    ncvar.attrib["cell_methods"] = "time: mean within years time: mean over years"

    #ncvar_deepest = def3D(ds,"$(varname)_deepest",ncvarattrib)
    #ncvar_deepest.attrib["long_name"] = "Deepest values of $(longname)"


    ncvar_Lx = Dict()
    ncvar_relerr = nothing

    if haskey(kw,:relerr)
        for (thresholds_name,thresholds_value) in thresholds

            ncvar_Lx[thresholds_value] = def4D(ds,"$(varname)_$(thresholds_name)",ncvarattrib)
            ncvar_Lx[thresholds_value].attrib["long_name"] = "$(longname) masked using relative error threshold $(thresholds_value)"

            #ncvar_deepest_Lx = def3D(ds,"$(varname)_deepest_$(thresholds_name)",ncvarattrib)
            #ncvar_deepest_Lx.attrib["long_name"] = "Deepest values of $(longname) masked using relative error threshold $(thresholds_value)"
        end

        # ncvar_err = def4D(ds,"$(varname)_err", type_save, ("lon", "lat", "depth", "time"))
        # ncvar_err.attrib["long_name"] = "Error standard deviation of $(longname)"
        # ncvar_err.attrib["units"] = units
        # ncvar_err.attrib["valid_min"] = type_save(0.0)
        # ncvar_err.attrib["valid_max"] = type_save(3.6)
        # ncvar_err.attrib["_FillValue"] = type_save(fillval)
        # ncvar_err.attrib["missing_value"] = type_save(fillval)

        ncvar_relerr = defVar(ds,"$(varname)_relerr", type_save, dimnames;
                              deflatelevel = deflatelevel,
                              chunksizes = chunksizes,
                              checksum = checksum,
                              )

        # section 3.1 'The conforming unit for quantities that represent fractions, or parts of a whole, is "1".'
        # https://web.archive.org/web/20171121154031/http://cfconventions.org/Data/cf-conventions/cf-conventions-1.6/build/cf-conventions.html
        ncvar_relerr.attrib["units"] = "1"
        ncvar_relerr.attrib["long_name"] = "Relative error of $(longname)"
        ncvar_relerr.attrib["valid_min"] = type_save(0.0)
        ncvar_relerr.attrib["valid_max"] = type_save(1.0)
        ncvar_relerr.attrib["_FillValue"] = type_save(fillval)
        ncvar_relerr.attrib["missing_value"] = type_save(fillval)
    end


    #    ncclimatology_bounds.attrib["climatology_bounds"] = type_save[244.0, 3622.0]

    # ncdatabins = defVar(ds,"databins", type_save, ("lon", "lat", "depth", "time");
    #                     deflatelevel = deflatelevel,
    #                     chunksizes = chunksizes)
    # ncdatabins.attrib["long_name"] = "Logarithm10 of number of data in bins"
    # ncdatabins.attrib["valid_min"] = type_save(0.0)
    # ncdatabins.attrib["valid_max"] = type_save(1.1)
    # ncdatabins.attrib['FillValue"] = type_save(fillval)
    # ncdatabins.attrib["missing_value"] = type_save(fillval)

    # ncoutlbins = defVar(ds,"outlbins", type_save, ("lon", "lat", "depth", "time"))
    # ncoutlbins.attrib["long_name"] = "Logarithm10 of number of outliers data in bins"
    # ncoutlbins.attrib["valid_min"] = type_save(0.0)
    # ncoutlbins.attrib["valid_max"] = type_save(fillval)
    # ncoutlbins.attrib["_FillValue"] = type_save(fillval)
    # ncoutlbins.attrib["missing_value"] = type_save(fillval)

    # Global attributes

    ds.attrib["Conventions"] = "CF-1.6"
    ds.attrib["title"] = "DIVA 4D analysis of $(longname)"
    ds.attrib["file_name"] = filename
    ds.attrib["product_id"] =  string(uuid1())
    ds.attrib["date"] = Dates.format(now(),"yyyy-mm-ddTHH:MM:SS")

    for (k,v) in ncglobalattrib
        ds.attrib[k] = v
    end


    # Define variables

    nclon[:]   = xyi[1][saveindex[1]]
    nclat[:]   = xyi[2][saveindex[2]]

    if idepth != -1
        ds["depth"][:]  = xyi[idepth][saveindex[3]]
    end

    if itime != -1
        ds["time"][:]  = xyi[itime]
    end

    # ncCLfield[:] = ...
    # ncCORRLEN[:] = ...
    # ncSNR[:] = ...
    # ncVARBACK[:] = ...

    if itime != -1 && haskey(kw,:climatology_bounds)
        ds["climatology_bounds"][:] = kw[:climatology_bounds]
    end

    sync(ds)
    return ncvar, ncvar_relerr, ncvar_Lx
end


"""
    ncvar, ncvar_relerr, ncvar_Lx, fi, relerr, index)

White a slice of data in a NetCDF given by the index `index`. The variable
`relerr` can be nothing.
"""
function writeslice(ncvar, ncvar_relerr, ncvar_Lx, fi, relerr, index;
                    saveindex = ntuple(i -> :, ndims(fi))
                    )
    fillval = NC_FILL_FLOAT

    tmp = copy(fi)
    tmp[isnan.(fi)] .= fillval
    ncvar[index...] = tmp[saveindex...]

    if relerr != nothing

        for (thresholds_value,ncvar_L) in ncvar_Lx
            tmp = copy(fi)
            tmp[isnan.(fi) .| (relerr .> thresholds_value)] .= fillval
            ncvar_L[index...] = tmp[saveindex...]
        end

        tmp = copy(relerr)
        tmp[isnan.(relerr)] .= fillval
        ncvar_relerr[index...] = tmp[saveindex...]
    end

    # ncvar_deepest[:] = ...
    # ncvar_deepest_L1[:] = ...
    # ncvar_deepest_L2[:] = ...
    # ncvar_err[:] = ...
    # ncdatabins[:] = ...
    # ncoutlbins[:] = ...
end


"""
    save(filename,xyi,fi,varname;
                          ncvarattrib = Dict(), ncglobalattrib = Dict(), ...)

Save the result of the analysis in a NetCDF file .

# Input arguments

* `filename`: the name of the NetCDF file
*  `xyi`: tuple with n vectors. Every element in this tuple represents a coordinate
  of the final grid on which the observations are interpolated
*  `fi`: the analysed field
* `varname`: the name of the NetCDF variable

# Optional arguments:
  * `ncglobalattrib`: a dictionary with the global attributes
  * `ncvarattrib`: a dictionary with the variable attributes
  * `relerr`: relative error

"""
function save(filename,xyi::NTuple{N,AbstractVector},fi,varname;
              kwargs...) where N


    kw = Dict(kwargs)
    # write the whole array
    index = (:,)

    Dataset(filename,"c") do ds
        ncvar, ncvar_relerr, ncvar_Lx = ncfile(ds,filename,xyi,varname;
                                               kwargs...)
        writeslice(ncvar, ncvar_relerr, ncvar_Lx,
                   fi, get(kw,:relerr,nothing), index)
    end

    return nothing
end


"""
    DIVAnd.saveobs(filename,xy,ids;
                   type_save = Float32,
                   timeorigin = DateTime(1900,1,1,0,0,0),
                   used = trues(size(ids)),
                   )
Save the location and time of the observation in the NetCDF file `filename` and
their identifier `ids`. `xy` is a tuple with the vectors longitude, latitude,
depth and time (as a vector of `DateTime`).

# Optional arguments:
  * `type_save`: the type to save the data (default Float32). However, the time
     is always saved as `Float64`.
  * `timeorigin`: time origin for the time units attribute (default is
1900-01-01 00:00:00)
  * `used`: allows to subset the data to save only used variables in the netCDF
     file

"""
function saveobs(filename,xy,ids;
                 type_save = Float32,
                 timeorigin = DateTime(1900,1,1,0,0,0),
                 used = trues(size(ids)),
                 checksum = :fletcher32,
                 chunksize = 10_000,
                 deflatelevel = 9,
                 )
    x,y,z,t = xy
    # keep only used observations
    xy = [xy_element[used] for xy_element in xy]
    ids = ids[used]

    # chunksizes should not exceed the number of observations
    chunksize = min(chunksize,length(ids))

    idlen = maximum(length.(ids))
    obsids = fill('\0',(idlen,length(ids)))
    for i = 1:length(ids)
        obsids[1:length(ids[i]),i] = Vector{Char}(ids[i])
    end

    mode = (isfile(filename) ? "a" : "c")

    Dataset(filename,mode) do ds
        #@show length(ids),idlen

        ds.dim["observations"] = length(ids)
        ds.dim["idlen"] = idlen

        ncobslon = defVar(ds,"obslon", type_save, ("observations",),
                          checksum = checksum,
                          deflatelevel = deflatelevel,
                          chunksizes = [chunksize],
                          )
        ncobslon.attrib["units"] = "degrees_east"
        ncobslon.attrib["standard_name"] = "longitude"
        ncobslon.attrib["long_name"] = "longitude"

        ncobslat = defVar(ds,"obslat", type_save, ("observations",),
                          checksum = checksum,
                          deflatelevel = deflatelevel,
                          chunksizes = [chunksize])
        ncobslat.attrib["units"] = "degrees_north"
        ncobslat.attrib["standard_name"] = "latitude"
        ncobslat.attrib["long_name"] = "latitude"

        ncobstime = defVar(ds,"obstime", Float64, ("observations",),
                           checksum = checksum,
                           deflatelevel = deflatelevel,
                           chunksizes = [chunksize])
        ncobstime.attrib["units"] = "days since " *
            Dates.format(timeorigin,"yyyy-mm-dd HH:MM:SS")

        ncobstime.attrib["standard_name"] = "time"
        ncobstime.attrib["long_name"] = "time"

        ncobsdepth = defVar(ds,"obsdepth", type_save, ("observations",),
                            checksum = checksum,
                            deflatelevel = deflatelevel,
                            chunksizes = [chunksize])

        ncobsdepth.attrib["units"] = "meters"
        ncobsdepth.attrib["positive"] = "down"
        ncobsdepth.attrib["standard_name"] = "depth"
        ncobsdepth.attrib["long_name"] = "depth below sea level"

        ncobsid = defVar(ds,"obsid", Char, ("idlen", "observations"),
                         checksum = checksum,
                         deflatelevel = deflatelevel,
                         chunksizes = [idlen,chunksize])
        ncobsid.attrib["long_name"] = "observation identifier"
        ncobsid.attrib["coordinates"] = "obstime obsdepth obslat obslon"

        ncobslon[:] = xy[1]
        ncobslat[:] = xy[2]
        ncobsdepth[:] = xy[3]
        # convertion is done in NCDatasets
        #ncobstime[:] = Dates.value.(Dates.Millisecond.(xy[4] - timeorigin)) / (24*60*60*1000.)
        ncobstime[:] = xy[4]
        ncobsid[:] = obsids
    end

    return nothing
end


"""
    DIVAnd.saveobs(filename,varname,value,xy,ids;
                   type_save = Float32,
                   timeorigin = DateTime(1900,1,1,0,0,0),
                   used = trues(size(ids)),
                   chunksize = 10_000,
                   )

Save `value` and the location and time of the observation in the NetCDF file `filename` and
their identifier `ids`. `xy` is a tuple with the vectors longitude, latitude,
depth and time (as a vector of `DateTime`). The values will be saved in the
variable called `varname`.

# Optional arguments:
  * `type_save`: the type to save the data (default Float32). However, the time
     is always saved as `Float64`.
  * `timeorigin`: time origin for the time units attribute (default is
1900-01-01 00:00:00)
  * `used`: allows to subset the data to save only used variables in the netCDF
     file

"""
function saveobs(filename,varname,value,xy,ids;
                 type_save = Float32,
                 timeorigin = DateTime(1900,1,1,0,0,0),
                 used = trues(size(ids)),
                 checksum = :fletcher32,
                 chunksize = 10_000,
                 deflatelevel = 9,
                 )

    # chunksizes should not exceed the number of observations
    chunksize = min(chunksize,sum(used))

    saveobs(filename,xy,ids;
            type_save = type_save,
            timeorigin = timeorigin,
            used = used,
            chunksize = chunksize,
            )


    Dataset(filename,"a") do ds
        ncobs = defVar(ds,varname, type_save, ("observations",),
                       checksum = checksum,
                       deflatelevel = deflatelevel,
                       chunksizes = [chunksize])
        ncobs[:] = value[used]
    end

    return nothing
end


@deprecate DIVAnd_save DIVAnd_save2
