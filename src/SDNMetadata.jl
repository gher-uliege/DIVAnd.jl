
const pathname = joinpath(dirname(@__FILE__), "..")

const PROJECTS = Dict(
    "EMODNET-chemistry" => Dict(
        "name" => "EMODnet Chemistry",
        "URL" => "https://emodnet.ec.europa.eu/en/chemistry",
        "baseurl_visualization" => "https://emodnet.ec.europa.eu/geoviewer",
        "baseurl_wms" => "http://ec.oceanbrowser.net/emodnet/Python/web/wms",
        "baseurl_http" => "https://emodnet.ec.europa.eu/geoviewer",
        "baseurl_opendap" => "http://opendap.oceanbrowser.net/thredds/dodsC/data/emodnet-domains",
        "template" => joinpath(pathname, "templates", "emodnet-chemistry.mustache"),
    ),
    "SeaDataNet" => Dict(
        "name" => "SeaDataNet",
        "URL" => "http://www.seadatanet.org/",
        "baseurl_visualization" => "http://sdn.oceanbrowser.net/web-vis/",
        "baseurl_wms" => "http://sdn.oceanbrowser.net/web-vis/Python/web/wms",
        "baseurl_http" => "http://sdn.oceanbrowser.net/data/SeaDataNet-domains",
        "baseurl_opendap" => "http://opendap.oceanbrowser.net/thredds/dodsC/data/SeaDataNet-domains",
        "template" => joinpath(pathname, "templates", "seadatanet.mustache"),
    ),
    "SeaDataCloud" => Dict(
        "name" => "SeaDataCloud",
        "URL" => "http://www.seadatanet.org/",
        "baseurl_visualization" => "http://sdn.oceanbrowser.net/web-vis/",
        "baseurl_wms" => "http://sdn.oceanbrowser.net/web-vis/Python/web/wms",
        "baseurl_http" => "http://sdn.oceanbrowser.net/data/SeaDataCloud-domains",
        "baseurl_opendap" => "http://opendap.oceanbrowser.net/thredds/dodsC/data/SeaDataNet-domains",
        "template" => joinpath(pathname, "templates", "seadatanet.mustache"),
    ),
)

const OriginatorEDMO_URL = "https://emodnet-chemistry.maris.nl/download/export.zip"


const layersep = "*"

"""
    escape(x)

Escape characters except slash.
"""
function escape(x)
    return join(HTTP.escapeuri.(split(x,"/")),"/")
end

"""encode parameters as key-value separated by : and +"""
encodeWMSStyle(params) = join([k * ':' * string(v) for (k, v) in params], "+")


"""
    ncglobalattrib,ncvarattrib = SDNMetadata(metadata,fi)

Based on the information in the dictionary `metadata` and the analysed 4D field
`fi` produce a list of NetCDF global and variable attributes for `DIVAnd_save2`.
To list all registered projects call `keys(DIVAnd.PROJECTS)`. For example:

```julia-repl
julia> using DIVAnd
julia> keys(DIVAnd.PROJECTS)
Base.KeySet for a Dict{String,Dict{String,String}} with 3 entries. Keys:
  "EMODNET-chemistry"
  "SeaDataNet"
  "SeaDataCloud"
```

"""
function SDNMetadata(
    metadata,
    filename,
    varname,
    lonr,
    latr;
    field = nothing,
    default_field_min = nothing,
    default_field_max = nothing,
    url_path = nothing,
)

    pathname = joinpath(dirname(@__FILE__), "..")

    sdn_parameter_urn = metadata["parameter_keyword_urn"]
    sdn_parameter = Vocab.resolve(sdn_parameter_urn)
    sdn_parameter_name = Vocab.prefLabel(sdn_parameter)

    # units
    sdn_uom = Vocab.findfirst(sdn_parameter, "related", "P06")
    sdn_uom_urn = Vocab.urn(sdn_uom)
    sdn_uom_name = Vocab.prefLabel(sdn_uom)


    # derived attributes

    project = PROJECTS[metadata["project"]]

    ncglobalattrib = OrderedDict{String,String}()

    institution = Vocab.resolve(metadata["institution_urn"])
    area = Vocab.resolve(metadata["area_keywords_urn"][1])
    parameter = Vocab.resolve(metadata["parameter_keyword_urn"])
    product_version = metadata["product_version"]
    product_type = "ANA"

    for (k, v) in metadata
        if k in ["netcdf_standard_name", "netcdf_long_name", "netcdf_units"]
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
            ncglobalattrib["search_keywords"] =
                join(Vocab.prefLabel.(Vocab.resolve.(v)), "|")
        elseif k == "area_keywords_urn"
            # Preferred label from SeaDataNet Vocabulary C19
            # example: Black Sea
            ncglobalattrib["area_keywords"] = join(Vocab.prefLabel.(Vocab.resolve.(v)), "|")
        elseif k == "product_version"
            # External shortname
            # example: SISMER-Atlantic Sea-Water_body_silicate-1.0-ANA
            ncglobalattrib["product_code"] = join(
                [
                    Vocab.abbreviated_name(institution),
                    Vocab.prefLabel(area),
                    Vocab.prefLabel(parameter),
                    product_version,
                    product_type,
                ],
                "-",
            )
        elseif k == "project"
            v = project["name"]
        end

        if typeof(v) <: Vector
            ncglobalattrib[k] = join(v, ", ")
        else
            ncglobalattrib[k] = v
        end

    end

    # Where the product can be downloaded
    ncglobalattrib["data_access"] = project["baseurl_http"]

    # Where the product can be visualised
    ncglobalattrib["WEB_visualisation"] = project["baseurl_visualization"]

    if (field != nothing) ||
       ((default_field_min != nothing) && (default_field_max != nothing))


        if (default_field_min == nothing) || (default_field_max == nothing)
            # default layer (first time instance and surface)
            default_field = field[:, :, end, 1]
            default_field_min, default_field_max =
                extrema(default_field[.!isnan.(default_field)])
        end

        domain = Vocab.prefLabel(area)
        preview_url = previewURL(
            filename,
            varname,
            lonr,
            latr,
            metadata["project"],
            domain;
            default_field_min = default_field_min,
            default_field_max = default_field_max,
            url_path = url_path,
        )


        # Example  http://ec.oceanbrowser.net/emodnet/Python/web/wms?styles=vmax%3A1.97872%2Bvmin%3A1.67568&format=image%2Fpng&height=500&bbox=27.5%2C42.0%2C30.0%2C44.5&decorated=true&transparent=true&layers=Back+Sea%2Fvarname.4Danl_autumn.nc%2Avarname&crs=CRS%3A84&service=WMS&request=GetMap&width=800&version=1.3.0
        ncglobalattrib["preview"] = preview_url
        @debug preview_url
    end

    ncvarattrib = OrderedDict(
        "units" => metadata["netcdf_units"],
        "standard_name" => metadata["netcdf_standard_name"],
        "long_name" => metadata["netcdf_long_name"],
        #  "sdn_parameter_urn" => sdn_parameter_urn,
        #  "sdn_parameter_name" => sdn_parameter_name,
        #  "sdn_uom_urn" => sdn_uom_urn,
        #  "sdn_uom_name" => sdn_uom_name,
    )

    return ncglobalattrib, ncvarattrib
end

"""
    contact = DIVAnd.getedmoinfo(edmo_code,role)

Returns a dictionary with the contact information from the EDMO registry
based on the prodivided `emdo_code`. `role` is the Sextant contact information
role, i.e. either "originator" or "author".
"""
function getedmoinfo(edmo_code, role)
    entry = DIVAnd.Vocab.EDMO()[edmo_code]

    contact = Dict{String,String}(
        "EDMO_CODE" => string(edmo_code),
        "EDMO_URL" => "https://edmo.seadatanet.org/report/$(edmo_code)",
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
    if endswith(fname, "zip")
        zp = ZipFile.Reader(fname)
        csvfile = zp.files[1]

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
    data::Array{String,2}, headers::Array{String,2} =
        readdlm(csvfile, ',', String; header = true)

    @debug "headers: $headers"
    @assert headers[:] == ["active", "author_edmo", "cdi_identifier", "originator_edmo"]


    # mapping from (author_edmo,cdi_identifier) to (active,originator_edmos)
    db = Dict{Tuple{Int64,String},Tuple{Bool,Vector{Int64}}}()

    for i = 1:size(data, 1)
        active, author_edmo, cdi_identifier, originator_edmos = data[i, :]
        key = (parse(Int, author_edmo), String(cdi_identifier))
        db[key] = (lowercase(active) == "true", parse.(Int64, split(originator_edmos, ',')))
    end

    return db
end


function writeerrors(notfound, errname)
    @info "Write error message in file: $(errname)"

    open(errname, "w") do f
        #log(0,"error message file opened")
        writedlm(f, reshape(["edmo", "local_cdi", "message"], 1, 3))

        # keep only unique rows
        data =
            sort(collect(Set([(v["edmo"], v["local_cdi"], v["status"]) for v in notfound])))
        #log(0,"write data")
        writedlm(f, data)
    end
end


function loadandcheckobsid(filepath)
    obsid = loadobsid(filepath)
    nbwrong = 0
    for i = 1:length(obsid)
        if !occursin('-', obsid[i])
            if nbwrong == 0
                @warn "in file $filepath observation identifier $(obsid[i]) does not contain a hypthen at index $i"
            end
            nbwrong += 1
        end
    end
    if nbwrong != 0
        error("$nbwrong invalid observation identifiers in $filepath")
    end
    return obsid
end

# EDMO information from originators

function getoriginators(
    db,
    filepaths::Vector{<:AbstractString},
    errname;
    ignore_errors = false,
)

    obsids = Iterators.flatten(loadandcheckobsid.(filepaths))

    originators, notfound =
        get_originators_from_obsid(db, obsids; ignore_errors = ignore_errors)

    if length(notfound) > 0
        writeerrors(notfound, errname)

        if ignore_errors
            @warn "Not all CDI could be found. See the file $(errname)"
        else
            error("Not all CDI could be found. See the file $(errname) ")
        end
    end

    return originators
end

function get_originators_from_obsid(db, obsids; ignore_errors = false)
    originators_edmo = Set{Int}()
    inactive = Dict{String,Any}[]
    notfound = Dict{String,Any}[]

    @info "Lookup obsids"

    for obsid in obsids
        (author_edmo_str, local_cdi) = split(obsid, '-'; limit = 2)
        author_edmo = parse(Int64, author_edmo_str)

        # strip trailing /v0 or /v1
        local_cdi = replace(local_cdi,r"/v\d+$" => "")

        key = (author_edmo, local_cdi)

        if haskey(db, key)
            v = db[key]
        else
            push!(
                notfound,
                Dict{String,Any}(
                    "edmo" => author_edmo,
                    "local_cdi" => local_cdi,
                    "status" =>
                        "CDI identifier is not present in CDI list ($(OriginatorEDMO_URL)). You might need to reload it.",
                ),
            )
            continue
        end

        #@show obsid,v
        active, originator_edmos_rec = v

        for oe in originator_edmos_rec
            push!(originators_edmo, oe)
        end

        if !active
            #@show inactive
            push!(
                inactive,
                Dict{String,Any}(
                    "edmo" => author_edmo,
                    "local_cdi" => local_cdi,
                    "status" => "inactive CDI",
                ),
            )
        end
    end

    @info "Query EDMO database"
    originators = Dict{String,String}[]
    for ae in collect(originators_edmo)
        originator = getedmoinfo(ae, "originator")

        if originator["name"] == "UNKNOWN"
            @info "ignoring originator: $(originator["name"])"
        else
            @info "originator: $(originator["name"])"
            push!(originators, originator)
        end
    end

    # sort by name
    sort!(originators, by = c -> c["name"])

    return originators, notfound
end

function labelandURL(s)
    c = DIVAnd.Vocab.resolve(s)
    return Dict("label" => DIVAnd.Vocab.prefLabel(c), "URL" => DIVAnd.Vocab.URL(c))
end

function URLsfromlabels(pname, labels)
    collection = DIVAnd.Vocab.SDNCollection(pname)
    concepts = Vocab.findbylabel(collection, labels)

    d = Dict{String,String}[]
    for i = 1:length(labels)
        c = concepts[i]
        push!(d, Dict("label" => DIVAnd.Vocab.prefLabel(c), "URL" => DIVAnd.Vocab.URL(c)))
    end

    return d
end

"""
    s = numstring(x)

Converts `x` to string and avoid to return "-0.0" if `x` is equal to zero.
https://github.com/gher-ulg/DIVAnd.jl/issues/28
"""
function numstring(x)
    s = string(x)
    if (s == "-0.0") || (x == 0)
        return "0"
    else
        return s
    end
end

function WMStimeparam(nctime, n)
    # heuristics as in
    # timemap as defined in divaonweb/divaonweb/util.py
    time = nctime[:]
    #time = [DateTime(2000,1,1),DateTime(2000,1,2)]; n = 1

    dtime = Dates.value(time[2] - time[1]) / (24*60*60*1000)
    seasonnames = ["winter","spring","summer","autumn"]
    month = Dates.month(nctime[n])
    year = Dates.year(nctime[n])
    season = (month-1) ÷ 3 + 1

    # heuristics
    if 27 <= dtime <= 32
        hr_t = @sprintf("%02d",month)
    elseif (27*3 <= dtime <= 32*3) && (length(time) <= 4)
        hr_t = seasonnames[n]
    elseif (27*3 <= dtime <= 32*3)
        hr_t = seasonnames[season] * " " * string(year)
    elseif 27*12 <= dtime <= 32*12
        hr_t = string(year)
    else
        hr_t = string(time[n])
    end

    @debug "WMS TIME parameter: $(string(nctime[n])); time increment $dtime; WMS parameter: $hr_t"
    return hr_t
end


function previewURL(
    filepath,
    varname_masked,
    lonr,
    latr,
    project,
    domain;
    colorbar_quantiles = [0.001, 0.999],
    basemap = "shadedrelief",
    preview_url_time = nothing,
    default_field_min = nothing,
    default_field_max = nothing,
    url_path = nothing,
)
    if url_path == nothing
        url_path = domain
    end

    preview_url_query = OrderedDict(
        "service" => "WMS",
        "request" => "GetMap",
        "version" => "1.3.0",
        "crs" => "CRS:84",
        "bbox" => "$(lonr[1]),$(latr[1]),$(lonr[end]),$(latr[end])",
        "decorated" => "true",
        "format" => "image/png",
        "layers" => url_path * "/" * filepath * layersep * varname_masked,
        "styles" => encodeWMSStyle(Dict(
            "vmin" => default_field_min,
            "vmax" => default_field_max,
        )),
        #"elevation" => "-0.0",
        #"time" => "03",
        "transparent" => "true",
        "height" => "500",
        "width" => "800",
        "basemap" => basemap,
    )

    if (preview_url_time != nothing)
        preview_url_query["time"] = preview_url_time
    end

    preview_url_query_string = string(HTTP.URI(; query = preview_url_query))

    # work-around for https://github.com/JuliaWeb/HTTP.jl/issues/323
    preview_url_query_string = replace(preview_url_query_string, r"^:" => "")

    preview_url = PROJECTS[project]["baseurl_wms"] * preview_url_query_string
    return preview_url
end


"""
    URL = previewURL(filepath,varname,project,domain;
                    colorbar_quantiles = [0.01,0.99],
                    basemap = "shadedrelief",
                    default_field_min = nothing,
                    default_field_max = nothing)

Generates a WMS preview URL
"""
function previewURL(
    filepath,
    varname,
    project,
    domain;
    colorbar_quantiles = [0.01, 0.99],
    basemap = "shadedrelief",
    default_field_min = nothing,
    default_field_max = nothing,
    url_path = nothing,
    default_error_level = "L2",
)

    Dataset(filepath, "r") do ds
        lon = ds["lon"][:]
        lonr = [minimum(lon), maximum(lon)]

        lat = ds["lat"][:]
        latr = [minimum(lat), maximum(lat)]

        # prefer layer L1 for plotting if available
        varname_masked = if haskey(ds, varname * "_" * default_error_level)
            varname * "_" * default_error_level
        else
            varname
        end

        # load default layer (first time instance and surface)
        # (time depth) lat lon
        var = ds[varname_masked]
        n_time = 1
        k_depth = 1
        min_non_missing = 0
        preview_url_time = nothing

        if ("time" in dimnames(var)) && ("depth" in dimnames(var))
            if size(var, 3) == 1
                # only one depth level
                field = var[:, :, 1, 1]
            else
                # find the time slice with the maximum data
                count_valid = zeros(Int, size(var, 4))
                for n = 1:size(var, 4)
                    field = var[:, :, k_depth, n]
                    count_valid[n] = sum(.!ismissing.(field))
                end
                @debug "count $(count_valid)"

                n_time = findmax(count_valid)[2]
                field = var[:, :, k_depth, n_time]
                if n_time != 1
                    preview_url_time = WMStimeparam(ds["time"], n_time)
                end
                #@debug "use n_time: $n_time $preview_url_time valid: $(count_valid[n_time])"
                @debug "n_time: $n_time"
                @debug "k_depth: $k_depth"
            end
        elseif "time" in dimnames(var)
            # only time
            field = var[:, :, 1]
        else
            # only depth
            if size(var, 3) == 1
                field = var[:, :, k_depth]
            else
                field = var[:, :, end]
            end
        end

        @debug "n_time: $n_time $(sum(.!ismissing.(field)))"
        if sum(.!ismissing.(field)) == 0
            @warn "only missing data in default view"
        end

        ranges = quantile(
            nomissing(field[.!(ismissing.(field) .| isnan.(field))]),
            colorbar_quantiles,
        )

        if default_field_min == nothing
            default_field_min = ranges[1]
        end

        if default_field_max == nothing
            default_field_max = ranges[2]
        end

        @debug "default_field_min $default_field_min"
        @debug "default_field_max $default_field_max"

        preview_url = previewURL(
            filepath,
            varname_masked,
            lonr,
            latr,
            project,
            domain;
            colorbar_quantiles = colorbar_quantiles,
            basemap = basemap,
            preview_url_time = preview_url_time,
            default_field_min = default_field_min,
            default_field_max = default_field_max,
            url_path = url_path,
        )

        return preview_url
    end
end

function gettemplatevars(
    filepaths::Vector{<:AbstractString},
    varname,
    project,
    cdilist;
    errname = split(filepaths[1], ".nc")[1] * ".cdi_import_errors.csv",
    WMSlayername = String[],
    previewindex = 1,
    basemap = "shadedrelief",
    additionalcontacts = [],
    ignore_errors = false,
    sigdigits = 5,
    url_path = nothing,
    WMSexclude = String[],
    default_error_level = "L2",
)

    # assume that grid and time coverage is the same as the
    # first file
    filepath = filepaths[1]
    isodateformat = DateFormat("yyyy-mm-ddTHH:MM:SS")

    baseurl_wms = PROJECTS[project]["baseurl_wms"]
    filename = basename(filepath)

    ds = Dataset(filepath, "r")
    lon = ds["lon"][:]
    lonr = [minimum(lon), maximum(lon)]
    dlon = (lonr[end] - lonr[1]) / (length(lon) - 1)

    lat = ds["lat"][:]
    latr = [minimum(lat), maximum(lat)]
    dlat = (latr[end] - latr[1]) / (length(lat) - 1)

    horizontal_resolution = if abs(dlon - dlat) < 1e-4
        dlat
    else
        @warn "warning: non uniform horizontal resolution $(dlon) $(dlat)"
        sqrt(dlon * dlat)
    end

    # round resolution
    horizontal_resolution = round(horizontal_resolution,sigdigits=sigdigits)

    # "degree (or km)",
    horizontal_resolution_units = "degree"

    date = ds.attrib["date"]
    date = replace(date, " " => "T")
    datetime = DateTime(date)

    product_id = ds.attrib["product_id"]

    # if multiple netCDF files are provided make sure that all netCDF
    # have the same product_id (global netCDF attribute)
    # as the first netCDF file
    for i = 2:length(filepaths)
        Dataset(filepaths[i], "a") do ds2
            ds2.attrib["product_id"] = product_id
        end
    end

    depth = if haskey(ds, "depth")
        ds["depth"][:]
    else
        [0.0]
    end

    doi = if haskey(ds.attrib, "doi")
        ds.attrib["doi"]
    else
        ""
    end

    # strip starting https://doi.org/ if necessary
    doi = replace(doi,r"^https://doi.org/" => "")

    temp_resolution_unit = "none"
    temp_resolution = "none"

    if haskey(ds, "time")
        nctime = ds["time"].var[:]

        dtime = if length(nctime) > 1
            nctime[2] - nctime[1]
        else
            zero(eltype(nctime))
        end

        timeunit = lowercase(ds["time"].attrib["units"])
        deltaunit, origin = split(timeunit, " since ")

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
            temp_resolution_unit = rstrip(deltaunit, 's')
            temp_resolution = nctime[2] - nctime[1]
        end
    end

    obstime = ds["obstime"][:]

    rmprefix(urn) = split(urn, ':')[end]

    P02_keywords = if haskey(ds.attrib, "search_keywords_urn")
        labelandURL.(split(ds.attrib["search_keywords_urn"], ", "))
    elseif haskey(ds.attrib, "search_keywords")
        URLsfromlabels("P02", split(ds.attrib["search_keywords"], '|'))
    else
        error("attribute search_keywords_urn or search_keywords are not found")
    end

    P35_keywords = if haskey(ds.attrib, "parameter_keyword_urn")
        labelandURL.(split(ds.attrib["parameter_keyword_urn"], ", "))
    elseif haskey(ds.attrib,"parameter_keyword")
        URLsfromlabels("P35", split(ds.attrib["parameter_keyword"], '|'))
    elseif haskey(ds.attrib,"parameter_keywords")
        URLsfromlabels("P35", split(ds.attrib["parameter_keywords"], '|'))
    else
        error("attribute parameter_keyword_urn or parameter_keyword(s) are not found")
    end

    C19_keywords = if haskey(ds.attrib, "area_keywords_urn")
        labelandURL.(split(ds.attrib["area_keywords_urn"], ", "))
    elseif haskey(ds.attrib, "area_keywords")
        URLsfromlabels("C19", split(ds.attrib["area_keywords"], '|'))
    else
        error("attribute area_keywords_urn or area_keywords are not found")
    end

    P36_urn = Vocab.urn.(Vocab.find(Vocab.Concept(P35_keywords[1]["URL"]), "broader", "P36"))
    P36_keywords = labelandURL.(P36_urn)

    @debug "P36_keywords: $(getindex.(P36_keywords,"label"))"

    domain = C19_keywords[1]["label"]

    edmo_code = if haskey(ds.attrib, "institution_urn")
        rmprefix.(ds.attrib["institution_urn"])
    elseif haskey(ds.attrib, "institution_edmo_code")
        ds.attrib["institution_edmo_code"]
    else
        error("attribute institution_urn or institution_edmo_code are not found")
    end

    product_code = get(ds.attrib, "product_code", "")


    # prefer layer L1/L2 for plotting if available
    varname_masked = if haskey(ds, varname * "_" * default_error_level)
        varname * "_" * default_error_level
    else
        varname
    end

    # random string for INSPIRE compliance
    time_period_id = randstring('a':'z', 32)

    templateVars = Dict(
        "project" => project,
        "product_id" => product_id,
        "product_code" => product_code,
        "product_version" => ds.attrib["product_version"],
        "abstract" => get(ds.attrib, "abstract", ""),
        "edmo_code" => edmo_code,
        "domain" => domain,
        "varname" => varname,
        "horizontal_resolution" => horizontal_resolution,
        "horizontal_resolution_units" => horizontal_resolution_units,
        "longitude_min" => minimum(lon),
        "longitude_max" => maximum(lon),
        "latitude_min" => minimum(lat),
        "latitude_max" => maximum(lat),
        "elevation_min" => numstring(minimum(-depth)),
        "elevation_max" => numstring(maximum(-depth)),
        "time_min" => Dates.format(minimum(obstime), isodateformat),
        "time_max" => Dates.format(maximum(obstime), isodateformat),
        "ndims" => ndims(ds[varname]),
        "nvertlevels" => length(depth),
        "temp_resolution_unit" => temp_resolution_unit,
        "temp_resolution" => temp_resolution,
        "CRS" => "WGS 84 (EPSG 4326)",
        "title" => ds.attrib["title"],
        "xml_creation_time" => Dates.format(now(), isodateformat),  # for XML Metadata
        "creation_date" => Dates.format(datetime, "yyyy-mm-dd"), # for the NetCDF file
        "update_date" => Dates.format(datetime, "yyyy-mm-dd"), # for the NetCDF file
        "netcdf_variables" => [], # added later
        "P02_keywords" => P02_keywords,
        "P35_keywords" => P35_keywords,
        "P36_keywords" => P36_keywords,
        "C19_keywords" => C19_keywords,
        "P02_date" => Dates.format(datetime, isodateformat),
        "P35_date" => Dates.format(datetime, isodateformat),
        "P36_date" => Dates.format(datetime, isodateformat),
        "C19_date" => Dates.format(datetime, isodateformat),
        "time_period_id" => time_period_id,
    )

    close(ds)
    #@show templateVars["elevation_min"]
    #@show templateVars["elevation_max"]

    for i = 1:length(filepaths)
        filepath_ = filepaths[i]
        Dataset(filepath_, "r") do ds
            for (name, var) in ds
                if (("lon" in dimnames(var)) && ("lat" in dimnames(var))) ||
                    (name == "obsid")

                    if name in WMSexclude
                        println("skipping ",name," for WMS layers")
                    else
                        if name == "obsid"
                            description = "Observations"
                        elseif haskey(var.attrib, "long_name")
                            description = var.attrib["long_name"]
                        else
                            description = name
                        end

                        # add WMS layer name suffix if provided
                        # this is useful if multiple NetCDF files are provided
                        if length(WMSlayername) >= i
                            if WMSlayername[i] != ""
                                description *= " ($(WMSlayername[i]))"
                            end
                        end

                        push!(templateVars["netcdf_variables"], (name, description, filepath_))
                    end
                end
            end
        end
    end

    DOI_URL = if doi != ""
        "https://doi.org/" * doi
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
    contacts = [getedmoinfo(parse(Int, edmo_code), "author")]

    db = loadoriginators(cdilist)

    originators = getoriginators(db, filepaths, errname, ignore_errors = ignore_errors)

    append!(contacts, originators)
    append!(contacts, additionalcontacts)

    bbox = join(
        [
            string(templateVars[c])
            for c in ["longitude_min", "latitude_min", "longitude_max", "latitude_max"]
        ],
        ',',
    )

    if url_path == nothing
        url_path = domain
    end

    preview_url =
        previewURL(filepaths[previewindex], varname, project, domain;
                   basemap = basemap,
                   url_path = url_path,
                   default_error_level = default_error_level,
                   )

    # Specify any input variables to the template as a dictionary.
    merge!(
        templateVars,
        Dict(
            "project_name" => PROJECTS[project]["name"],
            "project_url" => PROJECTS[project]["URL"],
            "preview" => preview_url,
            "L02_URL" => " http://vocab.nerc.ac.uk/collection/L02/current/006/",
            "L02_label" => "surface",
            "DOI_URL" => DOI_URL,
            # URL for the global data set
            "WMS_dataset_getcap" =>
                baseurl_wms * "?SERVICE=WMS&REQUEST=GetCapabilities&VERSION=1.3.0",
            "WMS_dataset_layer" => url_path * "/" * filepath * layersep * varname_masked,
            "WMS_layers" => [],
            "NetCDF_URL" => baseurl_http * "/" * escape(url_path) * "/" * escape(filepath),
            "NetCDF_description" => "Link to download the following file: " * filename,
            #        "OPENDAP_URL" => baseurl_opendap * "/" * domain * "/" * filepath * ".html",
            #        "OPENDAP_description" => "OPENDAP web page about the dataset " * filename,
            "OPENDAP_URLs" => [],
            "contacts" => contacts,
        ),
    )

    for i = 1:length(filepaths)
        # opendap
        description = "OPeNDAP web page about the dataset " * filename
        if length(WMSlayername) >= i
            if WMSlayername[i] != ""
                description *= " ($(WMSlayername[i]))"
            end
        end
        push!(
            templateVars["OPENDAP_URLs"],
            Dict(
                "URL" => baseurl_opendap * "/" * escape(url_path) * "/" * escape(filepaths[i]) * ".html",
                "description" => description,
            ),
        )
    end

    for (name, description, filepath_) in templateVars["netcdf_variables"]
        if (name == "obsid") || (endswith(name, "_L1") && !endswith(name, "deepest_L1"))
            layer_name = url_path * "/" * filepath_ * layersep * name
            if name == "obsid"
                layer_name = "point:" * layer_name
            end

            push!(
                templateVars["WMS_layers"],
                Dict(
                    "getcap" =>
                        baseurl_wms * "?SERVICE=WMS&REQUEST=GetCapabilities&VERSION=1.3.0",
                    "name" => layer_name,
                    "description" => "WMS layer for " * description,
                ),
            )
        end
    end

    #    @debug println("templateVars",templateVars)
    #    @debug println("templateVars NetCDF_URL",templateVars["NetCDF_URL"])
    @debug "WMS_layers: $(templateVars["WMS_layers"])"

    @info "Please check the following URLs"

    @info "URL (image preview): $preview_url"

    @info "NetCDF URL: $(templateVars["NetCDF_URL"])"

    for OPENDAP_URL in templateVars["OPENDAP_URLs"]
        @info "OPENDAP URL: $(OPENDAP_URL["URL"])"
    end

    for WMS_layer in templateVars["WMS_layers"]
        @debug "WMS layer name: $(WMS_layer["name"])"
    end

    return templateVars
end


function rendertemplate(templatefile, templateVars, xmlfilename)
    @info "Process template"

    template = read(templatefile, String)

    open(xmlfilename, "w") do f
        print(f, Mustache.render(template, templateVars))
    end
end



"""
    DIVAnd.divadoxml(filepath,varname,project,cdilist,xmlfilename;
                     ignore_errors = false,
                     WMSlayername = [],
                     previewindex = 1,
                     additionalcontacts = [],
                     additionalvars = Dict{String,Any}(),
                     default_error_level = "L2",
)

Generate the XML metadata file `xmlfilename` from the NetCDF
file `filepath` (or list of files) with the  NetCDF variable `varname`.
Project is either "SeaDataNet", "EMODNET-chemistry" or "SeaDataCloud".
`cdilist` is the file from $(OriginatorEDMO_URL).

The XML file contains a list of the data the originators. divadoxml
will abort with an error if some combinations of EDMO code, local CDI ID are
not present in the `cdilist`. Such errors can be ignored if `ignore_errors` is
set to true. To understand why some EDMO code/local CDI ID could not be found, one can
the decompress file from $(OriginatorEDMO_URL) which contains a file called
`export.csv`. This file has the columns `author_edmo` and `cdi_identifier` which this
function uses to find the data originators (column `originator_edmo`).

Information can be overridden with the dictionary `additionalvars`. The keys should
corresponds to the template tags found the in `template` directory. Template
tags are the strings inside {{ and }}.

NetCDF_URL should be suppplied since it's a URL of a ZIP file which is usually not from OceanBrowser.

If `filepath` is a vector of file names, the argument `WMSlayername` can be provided to give
additional information to distinguish between the NetCDF files. The elements of the vector of string
will be appended to the description of the WMS layer.

The resulting XML file includes the file names (provided by `filepath`).
Do not change the file names after running this function, otherwise the
XML will still contain a reference to the old file names. If you must change the
file names please do so before running this script.

If the data is present in a subfolder (e.g. "Winter") later on the OceanBrowser
webserver, the `filepath` should also contain this subfolder (e.g.
"Winter/somefile.nc"). The local directories should mirror the directory
structure on OceanBrowser. Relative paths should be used, and if the Julia code isn't right above the NetCDF
files, use cd("<path>") before each setting the files parameter which use paths relative to this path.

If the products will go into a sub-directory (which is different from the domain
name), the part of the path should with provided with the parameter `url_path`.

The link for the NetCDF download for EMODNET-chemistry for example is the following:

```
PROJECTS["EMODNET-chemistry"]["baseurl_http"] *  "/" * url_path * "/" * filepath
```

where `url_path` defaults to the the domain name (from the C19 vocabulary) if
it is not provided. The script will print the URLs for verification.

`additionalcontacts` is a list of dictionaries with additional condact
information to be added in the XML file. Elements are typically create by the
function `DIVAnd.getedmoinfo`.

`WMSexclude` is a list of string with NetCDF variables not be included the XML
under the WMS layer section.

`default_error_level` is either L1 or L2 corresponding to the masked variables
using relative error threshold 0.3 or 0.5. The default is L2.

### Example

```julia
download("$(OriginatorEDMO_URL)","export.zip")
files = [
         "Winter (January-March) - 6-year running averages/Water_body_chlorophyll-a.4Danl.nc",
         "Spring (April-June) - 6-year running averages/Water_body_chlorophyll-a.4Danl.nc",
         "Summer (July-September) - 6-year running averages/Water_body_chlorophyll-a.4Danl.nc",
         "Autumn (October-December) - 6-year running averages/Water_body_chlorophyll-a.4Danl.nc"
         ];

additionalcontacts = [
    DIVAnd.getedmoinfo(1977,"originator"), # US NODC for World Ocean Database
    DIVAnd.getedmoinfo(4630,"originator"), # CORIOLIS for CORA
]

DIVAnd.divadoxml(files,"Water_body_chlorophyll-a","EMODNET-chemistry","export.zip","test.xml";
    ignore_errors = true,
    additionalvars = Dict("abstract" => "Here goes the abstract"),
    additionalcontacts = additionalcontacts,
    WMSlayername = ["winter","spring","summer","autumn"]
)
```


For this function the following global NetCDF attributes are mandatory:

* `product_id`: UUID identifier
* `parameter_keyword_urn` or `parameter_keyword` (or `parameter_keywords`): P35 code (e.g. "SDN:P35::EPC00007") or preferred label (e.g. "Water body phosphate") respectively. Note the singular form (`parameter_keyword`) is preferred because there should be only a single P35 parameter per NetCDF file, but singular and plural forms are valid.
* `search_keywords_urn` or `search_keywords`: P02 code or preferred label
* `area_keywords_urn` or `area_keywords`: C19 code or preferred label
* `institution_edmo_code`: EDMO code number
* `product_version`: version of the product


Adding a global attribute, can be done using the following:

```julia
ds = NCDatafile("DIVA_file.nc","a")
ds.attrib["attribute_name"] = "attribute value"
close(ds)
```


"""
function divadoxml(
    filepaths::Vector{<:AbstractString},
    varname,
    project,
    cdilist,
    xmlfilename;
    ignore_errors = false,
    previewindex = 1,
    basemap = "shadedrelief",
    additionalvars = Dict{String,Any}(),
    additionalcontacts = [],
    WMSlayername = String[],
    WMSexclude = String[],
    url_path = nothing,
    default_error_level = "L2",
)

    # template file we will use.
    templatefile = PROJECTS[project]["template"]

    templateVars = gettemplatevars(
        filepaths,
        varname,
        project,
        cdilist,
        WMSlayername = WMSlayername,
        previewindex = previewindex,
        basemap = basemap,
        additionalcontacts = additionalcontacts,
        ignore_errors = ignore_errors,
        url_path = url_path,
        WMSexclude = WMSexclude,
        default_error_level = default_error_level,
    )

    merge!(templateVars, additionalvars)

    #@debug println("templateVars NetCDF_URL",templateVars["NetCDF_URL"])
    #@debug println("templateVars NetCDF_description",templateVars["NetCDF_description"])

    rendertemplate(templatefile, templateVars, xmlfilename)
end

function divadoxml(
    filepath::AbstractString,
    varname,
    project,
    cdilist,
    xmlfilename;
    WMSlayername = "",
    kwargs...,
)
    divadoxml(
        [filepath],
        varname,
        project,
        cdilist,
        xmlfilename;
        WMSlayername = [WMSlayername],
        kwargs...,
    )
end

"""
    divadoxml_originators(cdilist,filepaths,xmlfilename;
        errname = split(filepaths[1], ".nc")[1] * ".cdi_import_errors.csv",
        ignore_errors = false)

Produces the originators/contact XML fragment `xmlfilename` from the NetCDF files
in `filepaths`. See `divadoxml` for an explication of the optinal arguments.

## Example

```julia
using Glob
download(DIVAnd.OriginatorEDMO_URL,"export.zip")
run(`unzip export.zip`) # if unzip is installed
cdilist = "export.csv";
filepaths = glob("*/*/*silicate*nc")
xmlfilename = "out.xml"
ignore_errors = true
DIVAnd.divadoxml_originators(cdilist,filepaths,xmlfilename;
    ignore_errors = ignore_errors)
```
"""
function divadoxml_originators(
    cdilist,
    filepaths,
    xmlfilename;
    errname = split(filepaths[1], ".nc")[1] * ".cdi_import_errors.csv",
    templatefile = joinpath(
        dirname(pathof(DIVAnd)),
        "..",
        "templates",
        "emodnet-chemistry-originators.mustache",
    ),
    ignore_errors = false,
)

    ds = Dataset(filepaths[1], "r")
    edmo_code = if haskey(ds.attrib, "institution_urn")
        rmprefix.(ds.attrib["institution_urn"])
    else
        ds.attrib["institution_edmo_code"]
    end
    close(ds)

    db = loadoriginators(cdilist)
    originators = getoriginators(db, filepaths, errname, ignore_errors = ignore_errors)

    contacts = [getedmoinfo(parse(Int, edmo_code), "author")]
    append!(contacts, originators)

    templateVars = Dict("contacts" => contacts)

    rendertemplate(templatefile, templateVars, xmlfilename)
end

"""
    SDNObsMetadata(id)

Return a link to the SeaDataNet metadata page of the observation with the
identifier `id` (a combination of the EDMO code and local CDI ID).
This works only in IJulia.
"""
function SDNObsMetadata(id)
    edmo, local_CDI_ID = split(id, '-')
    url = "https://cdi.seadatanet.org/report/edmo/$(edmo)/$(local_CDI_ID)"

    display(
        "text/html",
        """
        Open in a new window <a target="blank" href="$(url)" >$(id)</a>
        <iframe width="900" height="700" src="$(url)"</iframe>
""",
    )
end



function updateNCfile(file; update_area_keywords = true)
    Dataset(file, "a") do ds
        if update_area_keywords
            area_keywords = ds.attrib["area_keywords"]
            # replace comma by the pipe symbol and strips white-space
            area_keywords = join(strip.(split(area_keywords, ",")), "|")
            @show area_keywords
            ds.attrib["area_keywords"] = area_keywords
        end
    end
end
