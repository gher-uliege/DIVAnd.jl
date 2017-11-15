function divand_save(filename,mask,varname,fi)

    sz = size(mask)

    ds = Dataset(filename,"c")

    # Dimensions

    ds.dim["time"] = sz[4];
    ds.dim["depth"] = sz[3];
    ds.dim["lat"] = sz[2];
    ds.dim["lon"] = sz[1];

    @show filename
    ncvar = defVar(ds, varname, Float32, ("lon", "lat", "depth", "time"))
    ncvar[:,:,:,:] = fi
    close(ds)

    return nothing
end


"""

ncinfo

```julia
ncglobalattrib = Dict(
    # Name of the project
    "project" => "SeaDataCloud",

    # Name from EDMO registry from http://seadatanet.maris2.nl/edmo/
    "institution" => "...",

    # Code number from EDMO registry
    "institution_edmo_code" => "...",

    # Production group
    "production" => "Diva group. E-mails: a.barth@ulg.ac.be, swatelet@ulg.ac.be",

    # Where the product can be downloaded
    "data_access" => "http://sdn.oceanbrowser.net/data/",

    # Where the product can be visualised
    "WEB_visualisation" => "http://sdn.oceanbrowser.net/web-vis/",

    # Name and emails from authors (seprated by comma)
    "Author_e-mail" => "Your Name1 <name1@example.com>, Your Name2 <name2@example.com>",

    # Source of the observation
    "source" => "observational data from SeaDataNet/EMODNet Chemistry Data Network",

    # Additional comment
    "comment" => "..."

    # Preferred label from SeaDataNet Parameter Discovery Vocabulary P02
    # example: Ammonium and ammonia concentration parameters in water bodies
    "search_keywords" => "...",

    # Preferred label from SeaDataNet Vocabulary P35
    # example: Water body ammonium
    "parameter_keywords" => "...",

    # Preferred label from SeaDataNet Vocabulary C19
    # example: Black Sea
    "area_keywords" => "...",

    "product_version" => "1.0",

    # Abstract for the product
    "abstract" => "...",

    # Example  http://ec.oceanbrowser.net/emodnet/Python/web/wms?styles=vmax%3A1.97872%2Bvmin%3A1.67568&format=image%2Fpng&height=500&bbox=27.5%2C42.0%2C30.0%2C44.5&decorated=true&transparent=true&layers=Back+Sea%2Fvarname.4Danl_autumn.nc%2Avarname&crs=CRS%3A84&service=WMS&request=GetMap&width=800&version=1.3.0
    "preview" => "...",

    # This option provides a place to acknowledge various types of support for the
    # project that produced the data
    "acknowledgment" => "...",

    # Digital Object Identifier of the data product
    "doi" => "...")
```

NetCDF CF Standard names
http://cfconventions.org/standard-names.html

[SeaDataNet Parameter Discovery Vocabulary P02](http://seadatanet.maris2.nl/v_bodc_vocab_v2/search.asp?lib=p02)

EMODnet Chemistry aggregated parameter names (P35)
http://seadatanet.maris2.nl/v_bodc_vocab_v2/search.asp?lib=p35

SeaVoX salt and fresh water body gazetteer (C19)
http://seadatanet.maris2.nl/v_bodc_vocab_v2/search.asp?lib=C19

EDMO registry
http://seadatanet.maris2.nl/edmo/
"""



function divand_save2(filename,mask,xyi,fi,varname; ncinfo = Dict(), ncglobalattrib = Dict())
    function defnD(ds,varname,dims,ncinfo)
        ncvar = defVar(ds,varname, Float32, dims)

        for a in ["units","standard_name","valid_min","valid_max"]
            if haskey(ncinfo,a)
                ncvar.attrib[a] = ncinfo[a]
            end
        end

        ncvar.attrib["_FillValue"] = Float32(fillval);
        ncvar.attrib["missing_value"] = Float32(fillval);

        return ncvar
    end

    def4D(ds,varname,ncinfo) = defnD(ds,varname,("lon", "lat", "depth", "time"),ncinfo)
    def3D(ds,varname,ncinfo) = defnD(ds,varname,("lon", "lat", "time"),ncinfo)


    (xi,yi,zi,ti) = xyi

    units = get(ncinfo,"units","")
    longname = get(ncinfo,"long_name","")
    validmin = get(ncinfo,"valid_min","")
    validmax = get(ncinfo,"valid_max",0.)

    sz = size(mask)

    fillval = NC_FILL_FLOAT

    ds = Dataset(filename,"c")

    # Dimensions

    ds.dim["time"] = sz[4];
    ds.dim["depth"] = sz[3];
    ds.dim["lat"] = sz[2];
    ds.dim["lon"] = sz[1];
    #ds.dim["nv"] = 2;
    #ds.dim["observations"] = 23626;
    #ds.dim["idlen"] = 40;

    # Declare variables

    # ncCLfield = defVar(ds,"CLfield", Float32, ("lon", "lat", "depth", "time"))
    # ncCLfield.attrib["long_name"] = "Correlation length field";
    # ncCLfield.attrib["valid_min"] = Float32(0.0);
    # ncCLfield.attrib["valid_max"] = Float32(0.76);
    # ncCLfield.attrib["_FillValue"] = Float32(fillval);
    # ncCLfield.attrib["missing_value"] = Float32(fillval);

    # ncCORRLEN = defVar(ds,"CORRLEN", Float32, ("depth", "time"))
    # ncCORRLEN.attrib["long_name"] = "Correlation Length";
    # ncCORRLEN.attrib["units"] = "degrees_north";

    # ncSNR = defVar(ds,"SNR", Float32, ("depth", "time"))
    # ncSNR.attrib["long_name"] = "Signal to Noise";

    # ncVARBACK = defVar(ds,"VARBACK", Float32, ("depth", "time"))
    # ncVARBACK.attrib["long_name"] = "Background Field Variance";
    # ncVARBACK.attrib["units"] = "umol/l^2";


    ncvar = def4D(ds,varname,ncinfo)
    ncvar.attrib["long_name"] = "$(longname)"
    ncvar.attrib["cell_methods"] = "time: mean within years time: mean over years";

    ncvar_L1 = def4D(ds,"$(varname)_L1",ncinfo)
    ncvar_L1.attrib["long_name"] = "$(longname) masked using relative error threshold 0.3";

    ncvar_L2 = def4D(ds,"$(varname)_L2",ncinfo)
    ncvar_L2.attrib["long_name"] = "$(longname) masked using relative error threshold 0.5";

    ncvar_deepest = def3D(ds,"$(varname)_deepest",ncinfo)
    ncvar_deepest.attrib["long_name"] = "Deepest values of $(longname)";

    ncvar_deepest_L1 = def3D(ds,"$(varname)_deepest_L1",ncinfo)
    ncvar_deepest_L1.attrib["long_name"] = "Deepest values of $(longname) masked using relative error threshold 0.3";

    ncvar_deepest_L2 = def3D(ds,"$(varname)_deepest_L2", ncinfo)
    ncvar_deepest_L2.attrib["long_name"] = "Deepest values of $(longname) masked using relative error threshold 0.5";

    # ncvar_err = def4D(ds,"$(varname)_err", Float32, ("lon", "lat", "depth", "time"))
    # ncvar_err.attrib["long_name"] = "Error standard deviation of $(longname)";
    # ncvar_err.attrib["units"] = units;
    # ncvar_err.attrib["valid_min"] = Float32(0.0);
    # ncvar_err.attrib["valid_max"] = Float32(3.6);
    # ncvar_err.attrib["_FillValue"] = Float32(fillval);
    # ncvar_err.attrib["missing_value"] = Float32(fillval);

    ncvar_relerr = defVar(ds,"$(varname)_relerr", Float32, ("lon", "lat", "depth", "time"))
    ncvar_relerr.attrib["long_name"] = "Relative error of $(longname)";
    ncvar_relerr.attrib["valid_min"] = Float32(0.0);
    ncvar_relerr.attrib["valid_max"] = Float32(1.0);
    ncvar_relerr.attrib["_FillValue"] = Float32(fillval);
    ncvar_relerr.attrib["missing_value"] = Float32(fillval);

    ncclimatology_bounds = defVar(ds,"climatology_bounds", Float32, ("nv", "time"))
    ncclimatology_bounds.attrib["climatology_bounds"] = Float32[244.0, 3622.0];

#     ncdatabins = defVar(ds,"databins", Float32, ("lon", "lat", "depth", "time"))
#     ncdatabins.attrib["long_name"] = "Logarithm10 of number of data in bins"
# ncdatabins.attrib["valid_min"] = Float32(0.0);
# ncdatabins.attrib["valid_max"] = Float32(1.1);
# ncdatabins.attrib["_FillValue"] = Float32(fillval);
# ncdatabins.attrib["missing_value"] = Float32(fillval);

ncdepth = defVar(ds,"depth", Float32, ("depth",))
ncdepth.attrib["units"] = "meters";
ncdepth.attrib["positive"] = "down";

nclat = defVar(ds,"lat", Float32, ("lat",))
nclat.attrib["units"] = "degrees_north";

nclon = defVar(ds,"lon", Float32, ("lon",))
nclon.attrib["units"] = "degrees_east";

ncobsdepth = defVar(ds,"obsdepth", Float32, ("observations",))
ncobsdepth.attrib["units"] = "meters";
ncobsdepth.attrib["positive"] = "down";

ncobsid = defVar(ds,"obsid", Char, ("idlen", "observations"))
ncobsid.attrib["long_name"] = "observation identifier";
ncobsid.attrib["coordinates"] = "obstime obsdepth obslat obslon";

ncobslat = defVar(ds,"obslat", Float32, ("observations",))
ncobslat.attrib["units"] = "degrees_north";

ncobslon = defVar(ds,"obslon", Float32, ("observations",))
ncobslon.attrib["units"] = "degrees_east";

ncobstime = defVar(ds,"obstime", Float32, ("observations",))
ncobstime.attrib["units"] = "days since 1900-01-01 00:00:00";

# ncoutlbins = defVar(ds,"outlbins", Float32, ("lon", "lat", "depth", "time"))
# ncoutlbins.attrib["long_name"] = "Logarithm10 of number of outliers data in bins";
# ncoutlbins.attrib["valid_min"] = Float32(0.0);
# ncoutlbins.attrib["valid_max"] = Float32(fillval);
# ncoutlbins.attrib["_FillValue"] = Float32(fillval);
# ncoutlbins.attrib["missing_value"] = Float32(fillval);

nctime = defVar(ds,"time", Float32, ("time",))
nctime.attrib["units"] = "Days since 1980-01-01";
nctime.attrib["climatology"] = "climatology_bounds";

# Global attributes

ds.attrib["Conventions"] = "CF-1.0";
ds.attrib["date"] = Dates.format(now(),"yyyy-mm-ddTHH:MM:SS")
ds.attrib["title"] = "DIVA 4D analysis of $(longname)";
ds.attrib["file_name"] = filename
ds.attrib["product_id"] =  Base.Random.uuid1()

for (k,v) in ncglobalattrib
    ncglobalattrib[k] = v
end


# Define variables

# ncCLfield[:] = ...
# ncCORRLEN[:] = ...
# ncSNR[:] = ...
# ncVARBACK[:] = ...
ncvar[:] = fi
# ncvar_L1[:] = ...
# ncvar_L2[:] = ...
# ncvar_deepest[:] = ...
# ncvar_deepest_L1[:] = ...
# ncvar_deepest_L2[:] = ...
# ncvar_err[:] = ...
# ncvar_relerr[:] = ...
# ncclimatology_bounds[:] = ...
# ncdatabins[:] = ...
# ncdepth[:] = ...
# nclat[:] = ...
# nclon[:] = ...
# ncobsdepth[:] = ...
# ncobsid[:] = ...
# ncobslat[:] = ...
# ncobslon[:] = ...
# ncobstime[:] = ...
# ncoutlbins[:] = ...
# nctime[:] = ...

close(ds)

# dims = [NcDim("longitude",sz[1]),
#         NcDim("latitude",sz[2]),
#         NcDim("depth",sz[3]),
#         NcDim("time",sz[4])]

# @show filename
# nc = NetCDF.create(filename,NcVar(varname,dims))
# nc[varname][:,:,:,:] = fi
# NetCDF.close(nc)

return nothing
end
