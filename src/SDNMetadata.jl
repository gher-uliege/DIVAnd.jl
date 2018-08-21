
const pathname = joinpath(dirname(@__FILE__),"..")

const PROJECTS = Dict(
    "EMODNET-chemistry" =>  Dict(
        "baseurl_visualization" => "http://ec.oceanbrowser.net/emodnet/",
        "baseurl_wms" => "http://ec.oceanbrowser.net/emodnet/Python/web/wms",
        "baseurl_http" => "http://ec.oceanbrowser.net/data/emodnet-domains",
        "baseurl_opendap" =>  "http://ec.oceanbrowser.net:8081/data/emodnet-domains",
        "template" => joinpath(pathname,"templates","emodnet-chemistry.mustache"),
    ),
    "SeaDataNet" => Dict(
        "baseurl_visualization" => "http://sdn.oceanbrowser.net/web-vis/",
        "baseurl_wms" => "http://sdn.oceanbrowser.net/web-vis/Python/web/wms",
        "baseurl_http" => "http://sdn.oceanbrowser.net/data/SeaDataNet-domains",
        "baseurl_opendap" => "http://sdn.oceanbrowser.net:8081/data/SeaDataNet-domains",
        "template" => joinpath(pathname,"templates","seadatanet.mustache"),
    ),
    "SeaDataCloud" => Dict(
        "baseurl_visualization" => "http://sdn.oceanbrowser.net/web-vis/",
        "baseurl_wms" => "http://sdn.oceanbrowser.net/web-vis/Python/web/wms",
        "baseurl_http" => "http://sdn.oceanbrowser.net/data/SeaDataCloud-domains",
        "baseurl_opendap" => "http://sdn.oceanbrowser.net:8081/data/SeaDataCloud-domains",
        "template" => joinpath(pathname,"templates","seadatanet.mustache"),
    )
)

const OriginatorEDMO_URL = "http://emodnet-chemistry.maris2.nl/download/export.zip"


const layersep = "*"

"""encode parameters as key-value separated by : and +"""
encodeWMSStyle(params) = join([k * ':' * string(v) for (k,v) in  params ],"+")


"""
    ncglobalattrib,ncvarattrib = SDNMetadata(metadata,fi)

Based on the information in the dictionary `metadata` and the analysed 4D field
`fi` produce a list of NetCDF global and variable attributes for `DIVAnd_save2`.
"""
function SDNMetadata(metadata,filename,varname,lonr,latr;
                     field = nothing,
                     default_field_min = nothing,
                     default_field_max = nothing)

    pathname = joinpath(dirname(@__FILE__),"..")

    sdn_parameter_urn = metadata["parameter_keyword_urn"]
    sdn_parameter = Vocab.resolve(sdn_parameter_urn)
    sdn_parameter_name = Vocab.prefLabel(sdn_parameter)

    # units
    sdn_uom = Vocab.findfirst(sdn_parameter,"related","P06")
    sdn_uom_urn = Vocab.urn(sdn_uom)
    sdn_uom_name = Vocab.prefLabel(sdn_uom)


    # derived attributes

    project = PROJECTS[metadata["project"]]

    ncglobalattrib = OrderedDict{String,String}()

    institution = Vocab.resolve(metadata["institution_urn"])
    area = Vocab.resolve(metadata["area_keywords_urn"][1])
    parameter = Vocab.resolve(metadata["parameter_keyword_urn"])
    product_version =  metadata["product_version"]
    product_type = "ANA"

    for (k,v) in metadata
        if k in ["netcdf_standard_name","netcdf_long_name","netcdf_units"]
            # these attribute do not need to be saved as global attributes
            continue
        end

        if k == "institution_urn"
            # Name from EDMO registry from http://seadatanet.maris2.nl/edmo/
            ncglobalattrib["institution"] = Vocab.name(institution)
        elseif k == "parameter_keyword_urn"
            # Preferred label from SeaDataNet Vocabulary P35
            # example: Water body ammonium
            ncglobalattrib["parameter_keyword"] = Vocab.prefLabel(parameter)
        elseif k == "search_keywords_urn"
            # Preferred label from SeaDataNet Parameter Discovery Vocabulary P02
            # example: Ammonium and ammonia concentration parameters in water bodies
            ncglobalattrib["search_keywords"] = join(Vocab.prefLabel.(Vocab.resolve.(v)),", ")
        elseif k == "area_keywords_urn"
            # Preferred label from SeaDataNet Vocabulary C19
            # example: Black Sea
            ncglobalattrib["area_keywords"] = join(Vocab.prefLabel.(Vocab.resolve.(v)),", ")
        elseif k == "product_version"
            # External shortname
            # example: SISMER-Atlantic Sea-Water_body_silicate-1.0-ANA
            ncglobalattrib["product_code"] = join([
                Vocab.abbreviated_name(institution),
                Vocab.prefLabel(area),
                Vocab.prefLabel(parameter),
                product_version,
                product_type],"-")
        end

        if typeof(v) <: Vector
            ncglobalattrib[k] = join(v,", ")
        else
            ncglobalattrib[k] = v
        end

    end

    # Where the product can be downloaded
    ncglobalattrib["data_access"] = project["baseurl_http"]

    # Where the product can be visualised
    ncglobalattrib["WEB_visualisation"] = project["baseurl_visualization"]

    if (field != nothing) || (
        (default_field_min != nothing) && (default_field_max != nothing))

        if (default_field_min == nothing) || (default_field_max == nothing)
            # default layer (first time instance and surface)
            default_field = field[:,:,end,1]
            default_field_min,default_field_max = extrema(default_field[.!isnan.(default_field)])
        end


        ncglobalattrib["preview"] = project["baseurl_wms"] * string(
            HTTP.URI(;query=
                     OrderedDict(
                         "service" => "WMS",
                         "request" => "GetMap",
                         "version" => "1.3.0",
                         "crs" => "CRS:84",
                         "bbox" => "$(lonr[1]),$(latr[1]),$(lonr[end]),$(latr[end])",
                         "decorated" => "true",
                         "format" => "image/png",
                         "layers" => Vocab.prefLabel(area) * "/" * filename * layersep * varname,
                         "styles" => encodeWMSStyle(Dict("vmin" => default_field_min,
                                                         "vmax" => default_field_max)),
                         #"elevation" => "-0.0",
                         #"time" => "03",
                         "transparent" => "true",
                         "height" => "500",
                         "width" => "800")))

        # Example  http://ec.oceanbrowser.net/emodnet/Python/web/wms?styles=vmax%3A1.97872%2Bvmin%3A1.67568&format=image%2Fpng&height=500&bbox=27.5%2C42.0%2C30.0%2C44.5&decorated=true&transparent=true&layers=Back+Sea%2Fvarname.4Danl_autumn.nc%2Avarname&crs=CRS%3A84&service=WMS&request=GetMap&width=800&version=1.3.0

    end

    ncvarattrib = OrderedDict("units" => metadata["netcdf_units"],
                              "standard_name" => metadata["netcdf_standard_name"],
                              "long_name" => metadata["netcdf_long_name"],
                              #  "sdn_parameter_urn" => sdn_parameter_urn,
                              #  "sdn_parameter_name" => sdn_parameter_name,
                              #  "sdn_uom_urn" => sdn_uom_urn,
                              #  "sdn_uom_name" => sdn_uom_name,
                              )

    return ncglobalattrib,ncvarattrib
end


function getedmoinfo(edmo_code,role)
    entry = DIVAnd.Vocab.EDMO()[edmo_code]

    contact = Dict{String,String}(
        "EDMO_CODE" => string(edmo_code),
        "EDMO_URL" => "http://seadatanet.maris2.nl/v_edmo/print.asp?n_code=$(edmo_code)",
        "name" => DIVAnd.Vocab.name(entry),
        "phone" => DIVAnd.Vocab.phone(entry),
        "fax" => DIVAnd.Vocab.fax(entry),
        "address" => DIVAnd.Vocab.address(entry),
        "city" => DIVAnd.Vocab.city(entry),
        "zip" => DIVAnd.Vocab.zipcode(entry),
        "country" => DIVAnd.Vocab.country(entry),
        "mail" => DIVAnd.Vocab.email(entry),
        "website" => DIVAnd.Vocab.website(entry),
        "role" => role,
    )

    return contact
end

"""
    db = loadoriginators(fname)

Load the CDI list from the file `fname`
(zip with a csv file, or csv file directly).
"""
function loadoriginators(fname::AbstractString)
    if endswith(fname,"zip")
        zp = ZipFile.Reader(fname);
        csvfile = zp.files[1];

        #db = loadoriginators(csvfile)

        buf = IOBuffer(read(csvfile))
        db = loadoriginators(buf)

        close(zp)

        return db
    else
        open(fname) do csvfile
            return loadoriginators(csvfile)
        end
    end
end

function loadoriginators(csvfile::IO)
    #data,headers = readdlm(csvfile,',', String; header = true) :: Tuple{Array{String,2},Vector{String}}
    data::Array{String,2}, headers::Array{String,2} = readdlm(csvfile,',', String; header = true)

    @assert headers[:] == ["active", "author_edmo",
                           "cdi_identifier", "originator_edmo"]


    # mapping from (author_edmo,cdi_identifier) to (active,originator_edmos)
    db = Dict{Tuple{Int64,String},Tuple{Bool,Vector{Int64}}}()

    for i = 1:size(data,1)
        active, author_edmo, cdi_identifier,originator_edmos = data[i,:]
        key = (parse(Int,author_edmo), String(cdi_identifier))
        db[key] = (lowercase(active) == "true",
                   parse.(Int64,split(originator_edmos,',')))
    end

    return db
end


function writeerrors(notfound,errname)
    @info "Write error message in file: $(errname)"

    open(errname,"w") do f
        #log(0,"error message file opened")
        writedlm(f,reshape(["edmo","local_cdi","message"],1,3))

        # keep only unique rows
        data = sort(collect(Set([(v["edmo"],v["local_cdi"],v["status"]) for v in notfound])))
        #log(0,"write data")
        writedlm(f,data)
    end
end

# EDMO information from originators

function getoriginators(db,filepaths::Vector{<:AbstractString},
                        errname; ignore_errors = false)

    obsids = Iterators.flatten(loadobsid.(filepaths))
    errname = replace(filepaths[1],r".nc$" => "") * ".cdi_import_errors_test.csv"
    originators,notfound = get_originators_from_obsid(
        db,obsids; ignore_errors = ignore_errors)

    if length(notfound) > 0
        writeerrors(notfound,errname)

        if ignore_errors
            @warn "Not all CDI could be found. See the file $(errname)"
        else
            error("Not all CDI could be found. See the file $(errname) ")
        end
    end

    return originators
end

function get_originators_from_obsid(db,obsids; ignore_errors = false)
    originators_edmo = Set{Int}()
    inactive = Dict{String,Any}[]
    notfound = Dict{String,Any}[]

    @info "Lookup obsids"

    for obsid in obsids
        (author_edmo_str,local_cdi)  = split(obsid,'-'; limit=2)
        author_edmo = parse(Int64,author_edmo_str)

        key = (author_edmo,local_cdi)

        if haskey(db,key)
            v = db[key]
        else
            push!(notfound,Dict{String,Any}(
                "edmo" => author_edmo,
                "local_cdi" => local_cdi,
                "status" => "CDI identifier is not present in CDI list ($(OriginatorEDMO_URL)). You might need to reload it."))
            continue
        end

        #@show obsid,v
        active,originator_edmos_rec  = v

        for oe in originator_edmos_rec
            push!(originators_edmo,oe)
        end

        if !active
            #@show inactive
            push!(inactive,Dict{String,Any}(
                "edmo" => author_edmo,
                "local_cdi" => local_cdi,
                "status" => "inactive CDI"))
        end
    end

    @info "Query EDMO database"
    originators = Dict{String,String}[]
    for ae in collect(originators_edmo)
        originator = getedmoinfo(ae,"originator")
        @info "originator: $(originator["name"])"
        push!(originators,originator)
    end

    # sort by name
    sort!(originators, by = c -> c["name"])

    return originators,notfound
end

function labelandURL(s)
    c = DIVAnd.Vocab.resolve(s)
    return Dict("label" => DIVAnd.Vocab.prefLabel(c),
                "URL" => DIVAnd.Vocab.URL(c))
end

function URLsfromlabels(pname,labels)
    collection = DIVAnd.Vocab.SDNCollection(pname)
    concepts = Vocab.findbylabel(collection,labels)

    d = Dict{String,String}[]
    for i = 1:length(labels)
        c = concepts[i]
        push!(d,Dict("label" => DIVAnd.Vocab.prefLabel(c),
                "URL" => DIVAnd.Vocab.URL(c)))
    end

    return d
end


function gettemplatevars(filepaths::Vector{<:AbstractString},varname,project,cdilist;
                         errname = split(filepaths[1],".nc")[1] * ".cdi_import_errors.csv",
                         WMSlayername = String[],
                         ignore_errors = false)

    # assume that grid and time coverage is the same as the
    # first file
    filepath = filepaths[1]
    isodateformat = DateFormat("yyyy-mm-ddTHH:MM:SS")

    baseurl_wms = PROJECTS[project]["baseurl_wms"]
    filename = basename(filepath)

    ds = Dataset(filepath,"r")
    lon = ds["lon"][:]
    lonr = [minimum(lon), maximum(lon)]
    dlon = (lonr[end]-lonr[1]) / (length(lon) - 1)

    lat = ds["lat"][:]
    latr = [minimum(lat), maximum(lat)]
    dlat = (latr[end]-latr[1]) / (length(lat) - 1)

    horizontal_resolution =
        if abs(dlon - dlat) < 1e-4
            dlat
        else
            @warn "warning: non uniform horizontal resolution $(dlon) $(dlat)"
            sqrt(dlon * dlat)
        end

    # "degree (or km)",
    horizontal_resolution_units = "degree"

    date = ds.attrib["date"]
    date = replace(date," " => "T")
    datetime = DateTime(date)

    product_id = ds.attrib["product_id"]

    depth =
        if haskey(ds,"depth")
            ds["depth"][:]
        else
            [0.]
        end

    doi =
        if haskey(ds.attrib,"doi")
            ds.attrib["doi"]
        else
            ""
        end

    temp_resolution_unit = "none"
    temp_resolution = "none"

    if haskey(ds,"time")
        nctime = ds["time"].var[:]
        dtime = nctime[2]-nctime[1]

        timeunit = lowercase(ds["time"].attrib["units"])
        deltaunit,origin = split(timeunit," since ")

        if length(nctime) == 1
            temp_resolution_unit = "none"
            temp_resolution = "none"
        elseif deltaunit == "months"
            if abs(dtime - 1) < 1e-6
                temp_resolution_unit = "month"
                temp_resolution = 1
            elseif abs(dtime - 3) < 1e-6
                temp_resolution_unit = "season"
                temp_resolution = 1
            elseif abs(dtime - 12) < 1e-6
                temp_resolution_unit = "year"
                temp_resolution = 1
            else
                temp_resolution_unit = "month"
                temp_resolution = dtime
            end
        else
            # remove plural s
            temp_resolution_unit = rstrip(deltaunit,'s')
            temp_resolution = nctime[2]-nctime[1]
        end
    end

    obstime = ds["obstime"][:]

    # load default layer (first time instance and surface)
    # (time depth) lat lon
    var = ds[varname]
    if ("time" in dimnames(var)) && ("depth" in dimnames(var))
        if size(var,3) == 1
            # only one depth level
            field = var[:,:,1,1]
        else
            field = var[:,:,end,1]
        end
    elseif "time" in dimnames(var)
        # only time
        field = var[:,:,1]
    else
        # only depth
        if size(var,3) == 1
            field = var[:,:,1]
        else
            field = var[:,:,end]
        end
    end

    default_field_min,default_field_max = extrema(field[.!ismissing.(field)])

    rmprefix(urn) = split(urn,':')[end]
    #P02_keywords = rmprefix.(split(ds.attrib["search_keywords_urn"]))
    #P35_keywords = rmprefix.(split(ds.attrib["parameter_keyword_urn"]))
    #C19_keywords = rmprefix.(split(ds.attrib["area_keywords_urn"]))

    P02_keywords =
        if haskey(ds.attrib,"search_keywords_urn")
            labelandURL.(split(ds.attrib["search_keywords_urn"]))
        else
            URLsfromlabels("P02",split(ds.attrib["search_keywords"],'|'))
        end

    P35_keywords =
        if haskey(ds.attrib,"parameter_keyword_urn")
            labelandURL.(split(ds.attrib["parameter_keyword_urn"]))
        else
            URLsfromlabels("P35",split(ds.attrib["parameter_keywords"],'|'))
        end

    C19_keywords =
        if haskey(ds.attrib,"area_keywords_urn")
            labelandURL.(split(ds.attrib["area_keywords_urn"]))
        else
            URLsfromlabels("C19",split(ds.attrib["area_keywords"],'|'))
        end

    domain = C19_keywords[1]["label"]

    edmo_code =
        if haskey(ds.attrib,"institution_urn")
            rmprefix.(ds.attrib["institution_urn"])
        else
            ds.attrib["institution_edmo_code"]
        end

    product_code = get(ds.attrib,"product_code","")

    templateVars = Dict(
        "project" => project,
        "product_id" => product_id,
        "product_code" => product_code,
        "product_version" => ds.attrib["product_version"],
        "abstract" => get(ds.attrib,"abstract",""),
        "edmo_code" => edmo_code,
        "domain" => domain,
        "varname" => varname,
        "horizontal_resolution" => horizontal_resolution,
        "horizontal_resolution_units" => horizontal_resolution_units,
        "longitude_min" => minimum(lon),
        "longitude_max" => maximum(lon),
        "latitude_min" => minimum(lat),
        "latitude_max" => maximum(lat),
        "elevation_min" => minimum(-depth),
        "elevation_max" => maximum(-depth),
        "time_min" => Dates.format(minimum(obstime),isodateformat),
        "time_max" => Dates.format(maximum(obstime),isodateformat),
        "default_field_min" => default_field_min,
        "default_field_max" => default_field_max,
        "ndims" => ndims(ds[varname]),
        "nvertlevels" => length(depth),
        "temp_resolution_unit" => temp_resolution_unit,
        "temp_resolution" => temp_resolution,
        "CRS" => "WGS 84 (EPSG 4326)",
        "title" => ds.attrib["title"],
        "xml_creation_time" => Dates.format(now(),isodateformat),  # for XML Metadata
        "creation_date" => Dates.format(datetime,"yyyy-mm-dd"), # for the NetCDF file
        "update_date" => Dates.format(datetime,"yyyy-mm-dd"), # for the NetCDF file
        "netcdf_variables" => [], # added later
        "P02_keywords" => P02_keywords,
        "P35_keywords" => P35_keywords,
        "C19_keywords" => C19_keywords,
        "P02_date" => Dates.format(datetime,isodateformat),
        "P35_date" => Dates.format(datetime,isodateformat),
        "C19_date" => Dates.format(datetime,isodateformat),
    )

    close(ds)

    for i = 1:length(filepaths)
        filepath = filepaths[i]
        Dataset(filepath,"r") do ds
            for (name,var) in ds
                if ("lon" in dimnames(var))  &&  ("lat" in dimnames(var))
                    if haskey(var.attrib,"long_name")
                        description = var.attrib["long_name"]
                    else
                        description = name
                    end

                    # add WMS layer name suffix if provided
                    # this is useful if multiple NetCDF files are provided
                    if length(WMSlayername) >= i
                        if WMSlayername[i] != ""
                            description *= " ($(WMSlayername[i]))";
                        end
                    end

                    push!(templateVars["netcdf_variables"],(name,description,filepath))
                end
            end
        end
    end

    DOI_URL =
        if doi != ""
            "http://dx.doi.org/" * doi
        else
            "na"
        end

    baseurl_wms = PROJECTS[project]["baseurl_wms"]
    baseurl_http = PROJECTS[project]["baseurl_http"]
    baseurl_opendap = PROJECTS[project]["baseurl_opendap"]

    @info "Loading EDMO information"

    # update 2018-08-09
    # old originator -> new author
    # old resource provider -> new originator
    contacts = [getedmoinfo(parse(Int,edmo_code),"author")]

    db = loadoriginators(cdilist)

    #UPDATE!! fixme
    #filepath = "/home/abarth/workspace/divadoxml-gui/Water_body_ammonium.4Danl_autumn.nc"

    originators = getoriginators(
        db,filepaths,errname,
        ignore_errors = ignore_errors)

    append!(contacts,originators)

    bbox = join([string(templateVars[c]) for c in [
        "longitude_min","latitude_min",
        "longitude_max","latitude_max"]],',')


    preview_url = PROJECTS[project]["baseurl_wms"] * string(
        HTTP.URI(;query=
                 OrderedDict(
                     "service" => "WMS",
                     "request" => "GetMap",
                     "version" => "1.3.0",
                     "crs" => "CRS:84",
                     "bbox" => "$(lonr[1]),$(latr[1]),$(lonr[end]),$(latr[end])",
                     "decorated" => "true",
                     "format" => "image/png",
                     "layers" => domain * "/" * filename * layersep * varname,
                     "styles" => encodeWMSStyle(Dict("vmin" => default_field_min,
                                                     "vmax" => default_field_max)),
                     #"elevation" => "-0.0",
                     #"time" => "03",
                     "transparent" => "true",
                     "height" => "500",
                     "width" => "800")))



    # Specify any input variables to the template as a dictionary.
    merge!(templateVars, Dict(
        "preview" => preview_url,
        "L02_URL" => " http://vocab.nerc.ac.uk/collection/L02/current/006/",
        "L02_label" => "surface",
        "DOI_URL" => DOI_URL,
        # URL for the global data set
        "WMS_dataset_getcap" => baseurl_wms * "?SERVICE=WMS&amp;REQUEST=GetCapabilities&amp;VERSION=1.3.0",
        "WMS_dataset_layer" => domain * "/" * filename * layersep * varname,
        "WMS_layers"  => [],
        "NetCDF_URL" => baseurl_http * "/" * domain * "/" * filename,
        "NetCDF_URL_description" => "Link to download the following file: " * filename,
        "OPENDAP_URL" => baseurl_opendap * "/" * domain * "/" * filename * ".html",
        "OPENDAP_description" => "OPENDAP web page about the dataset " * filename,
        "contacts" => contacts,
    ))

    for (name, description, filepath) in templateVars["netcdf_variables"]
        push!(templateVars["WMS_layers"],Dict(
            "getcap" => baseurl_wms * "?SERVICE=WMS&amp;REQUEST=GetCapabilities&amp;VERSION=1.3.0",
            "name" => domain * "/" * filepath * layersep * name,
            "description" => "WMS layer for " * description)
              )
    end

    return templateVars
end


function rendertemplate(templatefile,templateVars,xmlfilename)
    @info "Process template"

    template = read(templatefile,String)

    open(xmlfilename,"w") do f
        print(f,Mustache.render(template,templateVars))
    end
end



"""
    DIVAnd.divadoxml(filepath,varname,project,cdilist,xmlfilename;
                     ignore_errors = false,
                     WMSlayername = [],
                     additionalvars = Dict{String,Any}())

Generate the XML metadata file `xmlfilename` from the NetCDF
file `filepath` (or list of files) with the  NetCDF variable `varname`.
Project is either "SeaDataNet", "EMODNET-chemistry" or "SeaDataCloud".
`cdilist` is the file from $(OriginatorEDMO_URL).

The XML file contains a list of the data the originators. divadoxml
will abort with an error if some combinations of EDMO code, local CDI ID are
not present in the `cdilist`. Such errors can be ignored if `ignore_errors` is
set to true.

Information can be overriden with the dictionary `additionalvars`. The keys should
corresponds to the template tags found the in `template` directory. Template
tags are the strings inside {{ and }}.

If `filepath` is a vector of file names, the argument `WMSlayername` can be provided to give
additional information to distinguish between the NetCDF files. The elements of the vector of string
will be appended to the description of the WMS layer.

The resulting XML file includes the file names (provided by `filepath`).
Do not change the file names after running this function, otherwise the
XML will still contain a reference to the old file names. If you must change the
file names please do so before running this script.

### Example

```julia
download("$(OriginatorEDMO_URL)","export.zip")
files = [
         "Winter (January-March) - 6-year running averages/Water_body_chlorophyll-a.4Danl.nc",
         "Spring (April-June) - 6-year running averages/Water_body_chlorophyll-a.4Danl.nc",
         "Summer (July-September) - 6-year running averages/Water_body_chlorophyll-a.4Danl.nc",
         "Autumn (October-December) - 6-year running averages/Water_body_chlorophyll-a.4Danl.nc"
         ];


DIVAnd.divadoxml(files,"Water_body_chlorophyll-a","EMODNET-chemistry","export.zip","test.xml";
    ignore_errors = true,
    additionalvars = Dict("abstract" => "Here goes the abstract"),
    WMSlayername = ["winter","spring","summer","autumn"]
)
```
"""
function divadoxml(filepaths::Vector{<:AbstractString},varname,project,cdilist,xmlfilename;
                   ignore_errors = false,
                   additionalvars = Dict{String,Any}(), WMSlayername = String[])

    # template file we will use.
    templatefile = PROJECTS[project]["template"]

    templateVars = gettemplatevars(
        filepaths,varname,project,cdilist,
        WMSlayername = WMSlayername,
        ignore_errors = ignore_errors)

    merge!(templateVars,additionalvars)

    rendertemplate(templatefile,templateVars,xmlfilename)
end

function divadoxml(filepath::AbstractString,varname,project,cdilist,xmlfilename; WMSlayername = "", kwargs...)
    divadoxml([filepath],varname,project,cdilist,xmlfilename; WMSlayername = [WMSlayername], kwargs...)
end



"""
    SDNObsMetadata(id)

Return a link to the SeaDataNet metadata page of the observation with the
identifier `id` (a combination of the EDMO code and local CDI ID).
This works only in IJulia.
"""
function SDNObsMetadata(id)
    edmo,local_CDI_ID = split(id,'-')
    url = "http://seadatanet.maris2.nl/v_cdi_v3/print_wfs.asp" * string(
        HTTP.URI(;query=
                 OrderedDict(
                     "popup" => "yes",
                     "edmo" => edmo,
                     "identifier" => local_CDI_ID)))

    display("text/html", """
        Open in a new window <a target="blank" href="$(url)" >$(id)</a>
        <iframe width="900" height="700" src="$(url)"</iframe>
""")
end
