### DDR3 12800

#TODO: check keywords and mustach inverted sections
using Mustache
using NCDatasets
using DataArrays
using ZipFile
import divand
using DataStructures

include("../src/SDNMetadata.jl")
const layersep = "*"


function getedmoinfo(edmo_code,role)
    entry = divand.Vocab.EDMO()[edmo_code]

    contact = Dict(
        "EDMO_URL" => "http://seadatanet.maris2.nl/v_edmo/print.asp?n_code=$(edmo_code)",
        "name" => divand.Vocab.name(entry),
        "phone" => divand.Vocab.phone(entry),
        "fax" => divand.Vocab.fax(entry),
        "address" => divand.Vocab.address(entry),
        "city" => divand.Vocab.city(entry),
        "zip" => divand.Vocab.zipcode(entry),
        "country" => divand.Vocab.country(entry),
        "mail" => divand.Vocab.email(entry),
        "website" => divand.Vocab.website(entry),
        "role" => role,
    )

    return contact
end


#function loadNC(filepath,varname,project)


filepath = "/home/abarth/src/Diva-Workshops/notebooks/Water_body_Salinity.4Danl.nc"
varname = "Salinity"
project = "SeaDataCloud"

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
        warn("warning: non uniform horizontal resolution $(dlon) $(dlat)")
        sqrt(dlon * dlat)
    end

# "degree (or km)",
horizontal_resolution_units = "degree"

date = ds.attrib["date"]
date = replace(date," ","T")


product_id = ds.attrib["product_id"]

depth =
    if haskey(ds,"depth")
        ds["depth"][:]
    else
        [0.]
    end

temp_resolution_unit = "none"
temp_resolution = "none"

if haskey(ds,"time")
    time = ds["time"].var[:]
    dtime = time[2]-time[1]

    timeunit = lowercase(ds["time"].attrib["units"])
    deltaunit,origin = split(timeunit," since ")

    if length(time) == 1
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
        temp_resolution = time[2]-time[1]
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



P02_keywords = rmprefix.(split(ds.attrib["search_keywords_urn"]))
P35_keywords = rmprefix.(split(ds.attrib["parameter_keyword_urn"]))
C19_keywords = rmprefix.(split(ds.attrib["area_keywords_urn"]))

area = divand.Vocab.resolve(split(ds.attrib["area_keywords_urn"])[1])
domain = divand.Vocab.prefLabel(area)

edmo_code = rmprefix.(ds.attrib["institution_urn"])

product_code = ds.attrib["product_code"]

templatevar = Dict(
    "project" => project,
    "product_id" => product_id,
    "product_code" => product_code,
    "product_version" => ds.attrib["product_version"],
    "update_date" => date,
    "abstract" => ds.attrib["abstract"],
    "edmo_code" => edmo_code,
    "domain" => domain,
    "varname" => varname,
    "horizontal_resolution" => horizontal_resolution,
    "horizontal_resolution_units" => horizontal_resolution_units,
    "longitude_min" => minimum(lon),
    "longitude_max" => maximum(lon),
    "latitude_min" => minimum(lat),
    "latitude_max" => maximum(lat),
    "elevation_min" => minimum(depth),
    "elevation_max" => maximum(depth),
    "time_min" => string(minimum(obstime)),
    "time_max" => string(maximum(obstime)),
    "default_field_min" => default_field_min,
    "default_field_max" => default_field_max,
    # fix me
    "creation_time" => date,
    "ndims" => ndims(ds[varname]),
    "nvertlevels" => length(depth),
    "temp_resolution_unit" => temp_resolution_unit,
    "temp_resolution" => temp_resolution,
    "CRS" => "WGS 84 (EPSG 4326)",
    "title" => ds.attrib["title"],
    "creation_date" => ds.attrib["date"],
    "netcdf_variables" => [],
    "P02_keywords" => P02_keywords,
    "P35_keywords" => P35_keywords,
    "C19_keywords" => C19_keywords,
    "P02_date" => date,
    "P35_date" => date,
    "C19_date" => date,
)



for (name,var) in ds

    #print(name)
    if ("lon" in dimnames(var))  &&  ("lat" in dimnames(var))
        #print(name,var.ncattrs())
        if haskey(var.attrib,"long_name")
            description = var.attrib["long_name"]
        else
            description = name
        end

        push!(templatevar["netcdf_variables"],(name,description))

        # {"getcap": baseurl_wms + "?SERVICE=WMS&amp;REQUEST=GetCapabilities&amp;VERSION=1.3.0",
        #  "name": domain + "/" + filename + layersep + name,
        #  "description": "WMS layer for " + description});

    end
    #print(var.long_name)
end
close(ds)
#    return templatevar




user = copy(templatevar)




#    """producte the XML:
#      user: is a dict for the template variable
#    """

filename = basename(filepath)
project = user["project"]
# remove both lead and trailing whitespace
domain = strip(user["domain"])
edmo_code = user["edmo_code"]
varname = user["varname"]
#    divapath = user["divapath"]
doi = get(user,"doi","")

DOI_URL =
    if doi != ""
        "http://dx.doi.org/" * doi
    else
        "na"
    end

baseurl_wms = PROJECTS[project]["baseurl_wms"]
baseurl_http = PROJECTS[project]["baseurl_http"]
baseurl_opendap = PROJECTS[project]["baseurl_opendap"]

info("Loading EDMO information")

# entry = divand.Vocab.EDMO()[parse(Int,edmo_code)]

# contacts = [    Dict(
#             "EDMO_URL" => "http://seadatanet.maris2.nl/v_edmo/print.asp?n_code=" * edmo_code,
#             "name" => divand.Vocab.name(entry),
#             "phone" => divand.Vocab.phone(entry),
#             "fax" => divand.Vocab.fax(entry),
#             "address" => divand.Vocab.address(entry),
#             "city" => divand.Vocab.city(entry),
#             "zip" => divand.Vocab.zipcode(entry),
#             "country" => divand.Vocab.country(entry),
#             "mail" => divand.Vocab.email(entry),
#             "website" => divand.Vocab.website(entry),
#             "role" => "originator",
#         )]

contacts = [getedmoinfo(parse(Int,edmo_code),"originator")]

# EDMO information from originators
errname = split(filepath,".nc")[1] * ".cdi_import_errors.csv"


fname = "/home/abarth/workspace/divadoxml-gui/CDI-list-export.zip"
zp = ZipFile.Reader(fname); csvfile = zp.files[1];
data,headers = readdlm(csvfile,',', String; header = true);

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


filepath = "/home/abarth/workspace/divadoxml-gui/Water_body_ammonium.4Danl_autumn.nc"
nc = Dataset(filepath,"r")
obsids = nc["obsid"][:]
close(nc)


originators_edmo = Set{Int}()
inactive = []
notfound = []

OriginatorEDMO_URL = "http://emodnet-chemistry.maris2.nl/download/export.zip"

errname = replace(filepath,r".nc$","") * ".cdi_import_errors_test.csv"
ignore_errors = true



info("Lookup obsids")

for i = 1:size(obsids,2)
    obsid = strip(join(obsids[:,i]),'\0')

    (author_edmo_str,local_cdi)  = split(obsid,'-'; limit=2)
    author_edmo = parse(Int64,author_edmo_str)

    key = (author_edmo,local_cdi)

    if haskey(db,key)
        v = db[key]
    else
        push!(notfound,
              Dict("edmo" => author_edmo, "local_cdi" => local_cdi,
                   "status" => "CDI identifier is not present in CDI list ($(OriginatorEDMO_URL)). You might need to reload it."))
        continue
    end

    #@show obsid,v
    active,originator_edmos_rec  = v

    for oe in originator_edmos_rec
        push!(originators_edmo,oe)
    end

    if !active
        append!(inactive,Dict("edmo" => author_edmo,
                              "local_cdi" => local_cdi,
                              "status" => "inactive CDI"))
    end
end

#log(0,"originators_edmo: %s" % (str(originators_edmo),))

if length(notfound) > 0
    #log(0,"write error message in file: %s" % (errname,))


    open(errname,"w") do f
        #log(0,"error message file opened")
        writedlm(f,reshape(["edmo","local_cdi","message"],1,3))

        # keep only unique rows
        data = sort(collect(Set([(v["edmo"],v["local_cdi"],v["status"]) for v in notfound])))
        #log(0,"write data")
        writedlm(f,data)
    end

    if ignore_errors
        warn("Not all CDI could be found. See the file $(errname)")
    else
        error("Not all CDI could be found. See the file $(errname) ")
    end

end

info("Query EDMO database")
originators = []
for ae in sort(collect(originators_edmo))
    originator = getedmoinfo(ae,"resourceProvider")
    println("resource provider: $(originator["name"])")
    push!(originators,originator)
end




#=
originators = get_originators_from_obsid(filepath,errname,user["ignore_errors"])
=#


append!(contacts,originators)

info("Process template")
# template file we will use.
templatefile = PROJECTS[project]["template"]


bbox = join([string(user[c]) for c in ["longitude_min","latitude_min","longitude_max","latitude_max"]],',')
#print("entry",entry,bbox,user["default_field_max"])

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


#print("preview",preview_url);


# Specify any input variables to the template as a dictionary.
templateVars = Dict(
    "preview" => preview_url,
    "L02_URL" => " http://vocab.nerc.ac.uk/collection/L02/current/006/",
    "L02_label" => "surface",
    "DOI_URL" => DOI_URL,
    # URL for the global data set
    "WMS_dataset_getcap" => baseurl_wms * "?SERVICE=WMS&amp;REQUEST=GetCapabilities&amp;VERSION=1.3.0",
    "WMS_dataset_layer" => domain * "/" * filename * layersep * varname,
    # URL for specific layers
    # "WMS_layers"  => [{"getcap" => baseurl_wms * "?SERVICE=WMS&amp;REQUEST=GetCapabilities&amp;VERSION=1.3.0",
    #                 "name" => domain * "/" * filename * layersep * varname,
    #                 "description" => "WMS layer for " * user["title"]}],
    "WMS_layers"  => [],
    "NetCDF_URL" => baseurl_http * "/" * domain * "/" * filename,
    "NetCDF_URL_description" => "Link to download the following file: " * filename,
    "OPENDAP_URL" => baseurl_opendap * "/" * domain * "/" * filename * ".html",
    "OPENDAP_description" => "OPENDAP web page about the dataset " * filename,
    "contacts" => contacts,
)


for (name, description) in user["netcdf_variables"]
    push!(templateVars["WMS_layers"],Dict(
        "getcap" => baseurl_wms * "?SERVICE=WMS&amp;REQUEST=GetCapabilities&amp;VERSION=1.3.0",
        "name" => domain * "/" * filename * layersep * name,
        "description" => "WMS layer for " * description))
end


merge!(templateVars,user)


templatefile = "/home/abarth/workspace/divadoxml-gui/templates/emodnet-chemistry.mustache"
template = readstring(open(templatefile))

xmlfilename = "test.xml"

open(xmlfilename,"w") do f
    print(f,Mustache.render(template,templateVars))
end

