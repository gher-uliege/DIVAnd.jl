function divand_save(filename,mask,varname,fi)

    sz = size(mask)

    ds = Dataset(filename,"c")

    # Dimensions

    ds.dim["time"] = sz[4]
    ds.dim["depth"] = sz[3]
    ds.dim["lat"] = sz[2]
    ds.dim["lon"] = sz[1]

    @show filename
    ncvar = defVar(ds, varname, Float32, ("lon", "lat", "depth", "time"))
    ncvar[:,:,:,:] = fi
    close(ds)

    return nothing
end


"""
    divand_save2(filename,mask,xyi,fi,varname;
                      ncvarattrib = Dict(), ncglobalattrib = Dict(), ...)

Save the result of the analysis in a NetCDF file .

# Input arguments

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

"""



function divand_save2(filename,mask,xyi,fi,varname;
                      ncvarattrib = Dict(), ncglobalattrib = Dict(),
                      thresholds = [("L1",0.3),("L2",0.5)],
                      deflatelevel = 5,
                      chunksizes = [100,100,1,1],
                      type_save = Float32,
                      kwargs...)

    function defnD(ds,varname,dims,ncvarattrib)
        ncvar = defVar(ds,varname, type_save, dims;
                       deflatelevel = deflatelevel,
                       chunksizes = chunksizes[1:length(dims)])

        for (k,v) in ncvarattrib
            ncvar.attrib[k] = ncvarattrib[v]
        end

        ncvar.attrib["_FillValue"] = type_save(fillval)
        ncvar.attrib["missing_value"] = type_save(fillval)

        return ncvar
    end

    def4D(ds,varname,ncvarattrib) = defnD(ds,varname,("lon", "lat", "depth", "time"),ncvarattrib)
    def3D(ds,varname,ncvarattrib) = defnD(ds,varname,("lon", "lat", "time"),ncvarattrib)

    kw = Dict(kwargs)
    # chunksizes should not exceed the size of fi
    chunksizes = min.(chunksizes,collect(size(fi)))
    @show chunksizes

    (xi,yi,zi,ti) = xyi

    units = get(ncvarattrib,"units","")
    longname = get(ncvarattrib,"long_name","")
    validmin = get(ncvarattrib,"valid_min","")
    validmax = get(ncvarattrib,"valid_max",0.)

    sz = size(mask)

    fillval = NC_FILL_FLOAT

    Dataset(filename,"c") do ds

        # Dimensions

        ds.dim["lon"] = sz[1]
        ds.dim["lat"] = sz[2]
        ds.dim["depth"] = sz[3]
        ds.dim["time"] = sz[4]
        ds.dim["nv"] = 2

        @show sz
        # Declare variables

        nclon = defVar(ds,"lon", Float64, ("lon",))
        nclon.attrib["units"] = "degrees_east"
        nclon.attrib["standard_name"] = "longitude"
        nclon.attrib["long_name"] = "longitude"


        nclat = defVar(ds,"lat", Float64, ("lat",))
        nclat.attrib["units"] = "degrees_north"
        nclat.attrib["standard_name"] = "latitude"
        nclat.attrib["long_name"] = "latitude"

        ncdepth = defVar(ds,"depth", Float64, ("depth",))
        ncdepth.attrib["units"] = "meters"
        ncdepth.attrib["positive"] = "down"
        ncdepth.attrib["standard_name"] = "depth"
        ncdepth.attrib["long_name"] = "depth below sea level"

        nctime = defVar(ds,"time", Float64, ("time",))
        nctime.attrib["units"] = "days since 1900-01-01 00:00:00"
        nctime.attrib["standard_name"] = "time"
        nctime.attrib["long_name"] = "time"
        nctime.attrib["calendar"] = "standard"

        if haskey(kw,:climatology_bounds)
            nctime.attrib["climatology"] = "climatology_bounds"
            ncclimatology_bounds = defVar(ds,"climatology_bounds", Float64, ("nv", "time"))
            ncclimatology_bounds.attrib["units"] = "days since 1900-01-01 00:00:00"
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

        ncvar_deepest = def3D(ds,"$(varname)_deepest",ncvarattrib)
        ncvar_deepest.attrib["long_name"] = "Deepest values of $(longname)"



        if haskey(kw,:relerr)
            for (thresholds_name,thresholds_value) in thresholds

                ncvar_Lx = def4D(ds,"$(varname)_$(thresholds_name)",ncvarattrib)
                ncvar_Lx.attrib["long_name"] = "$(longname) masked using relative error threshold $(thresholds_value)"

                ncvar_deepest_Lx = def3D(ds,"$(varname)_deepest_$(thresholds_name)",ncvarattrib)
                ncvar_deepest_Lx.attrib["long_name"] = "Deepest values of $(longname) masked using relative error threshold $(thresholds_value)"
            end

            # ncvar_err = def4D(ds,"$(varname)_err", type_save, ("lon", "lat", "depth", "time"))
            # ncvar_err.attrib["long_name"] = "Error standard deviation of $(longname)"
            # ncvar_err.attrib["units"] = units
            # ncvar_err.attrib["valid_min"] = type_save(0.0)
            # ncvar_err.attrib["valid_max"] = type_save(3.6)
            # ncvar_err.attrib["_FillValue"] = type_save(fillval)
            # ncvar_err.attrib["missing_value"] = type_save(fillval)

            @show "relerr"
            ncvar_relerr = defVar(ds,"$(varname)_relerr", type_save, ("lon", "lat", "depth", "time");
                                  deflatelevel = deflatelevel,
                                  chunksizes = chunksizes)

            @show "relerr2"
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
        ds.attrib["product_id"] =  repr(Base.Random.uuid1())
        ds.attrib["date"] = Dates.format(now(),"yyyy-mm-ddTHH:MM:SS")

        for (k,v) in ncglobalattrib
            ds.attrib[k] = v
        end


        # Define variables

        nclon[:]   = xyi[1]
        nclat[:]   = xyi[2]
ncdepth[:] = xyi[3]
nctime[:]  = xyi[4]

# ncCLfield[:] = ...
# ncCORRLEN[:] = ...
# ncSNR[:] = ...
# ncVARBACK[:] = ...
ncvar[:] = DataArray(fi,isnan.(fi))

if haskey(kw,:relerr)
    relerr = kw[:relerr]

    for (thresholds_name,thresholds_value) in thresholds
        ds["$(varname)_$(thresholds_name)"][:] =
            DataArray(fi,isnan.(fi) .| (relerr .> thresholds_value))
    end

    ncvar_relerr[:] = DataArray(relerr,isnan.(fi))
end

# ncvar_deepest[:] = ...
# ncvar_deepest_L1[:] = ...
# ncvar_deepest_L2[:] = ...
# ncvar_err[:] = ...
if haskey(kw,:climatology_bounds)
    ncclimatology_bounds[:] = kw[:climatology_bounds]
end
# ncdatabins[:] = ...
# ncoutlbins[:] = ...
end


return nothing
end


function divand_save_obs(filename,ids,xy; type_save = Float32)
    x,y,z,t = xy

    idlen = maximum(length.(ids))
    obsids = cat(2,[convert(Vector{Char},rpad(id,idlen)) for id in ids]...)

    Dataset(filename,"a") do

        ds.dim["observations"] = length(ids)
        ds.dim["idlen"] = idlen

        ncobslon = defVar(ds,"obslon", type_save, ("observations",))
        ncobslon.attrib["units"] = "degrees_east"
        ncobslon.attrib["standard_name"] = "longitude"
        ncobslon.attrib["long_name"] = "longitude"

        ncobslat = defVar(ds,"obslat", type_save, ("observations",))
        ncobslat.attrib["units"] = "degrees_north"
        ncobslat.attrib["standard_name"] = "latitude"
        ncobslat.attrib["long_name"] = "latitude"

        ncobstime = defVar(ds,"obstime", Float64, ("observations",))
        ncobstime.attrib["units"] = "days since 1900-01-01 00:00:00"
        ncobstime.attrib["standard_name"] = "time"
        ncobstime.attrib["long_name"] = "time"

        ncobsdepth = defVar(ds,"obsdepth", type_save, ("observations",))
        ncobsdepth.attrib["units"] = "meters"
        ncobsdepth.attrib["positive"] = "down"
        ncobsdepth.attrib["standard_name"] = "depth"
        ncobsdepth.attrib["long_name"] = "depth below sea level"

        ncobsid = defVar(ds,"obsid", Char, ("idlen", "observations"))
        ncobsid.attrib["long_name"] = "observation identifier"
        ncobsid.attrib["coordinates"] = "obstime obsdepth obslat obslon"

        ncobslon[:] = xy[1]
        ncobslat[:] = xy[2]
        ncobsdepth[:] = xy[3]
        ncobstime[:] = xy[4]
        ncobsid[:] = obsids
    end
end


@deprecate divand_save divand_save2
