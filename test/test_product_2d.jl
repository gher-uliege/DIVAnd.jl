import DIVAnd
if VERSION >= v"0.7.0-beta.0"
    using Test
else
    using Base.Test
end
using DataStructures
using Missings
using NCDatasets

#
varname = "Salinity"
filename = "WOD-Salinity.nc"


bathname = joinpath(dirname(@__FILE__),"..","..","DIVAnd-example-data",
                    "Global","Bathymetry","gebco_30sec_16.nc")
bathisglobal = true

obsname = joinpath(dirname(@__FILE__),"..","..","DIVAnd-example-data",
                    "Provencal","WOD-Salinity.nc")

cdilist = joinpath(dirname(@__FILE__),"..","data","CDI-list-export.csv")


if !isfile(bathname)
    bathname = download("https://dox.ulg.ac.be/index.php/s/U0pqyXhcQrXjEUX/download")
end

if !isfile(obsname)
    obsname = download("https://dox.ulg.ac.be/index.php/s/PztJfSEnc8Cr3XN/download")
end

obsname = joinpath(dirname(@__FILE__),"..","data","sample-file.nc")

obsvalue,obslon,obslat,obsdepth,obstime,obsids = DIVAnd.loadobs(Float64,obsname,"Salinity")

sel = obsdepth .< 10

dx = dy = 0.5
lonr = 3:dx:11.8
latr = 42.0:dy:44.0
epsilon2 = 0.01


sz = (length(lonr),length(latr))

lenx = fill(200_000,sz)
leny = fill(200_000,sz)

years = 1993:1993

year_window = 10

# winter: January-March    1,2,3
# spring: April-June       4,5,6
# summer: July-September   7,8,9
# autumn: October-December 10,11,12

monthlists = [
    [1,2,3],
    [4,5,6],
    [7,8,9],
    [10,11,12]
];


TS = DIVAnd.TimeSelectorYW(years,year_window,monthlists)

varname = "Salinity"

# File name based on the variable (but all spaces are replaced by _)
filename = "Water_body_$(replace(varname,' ' => '_')).3Danl.nc"

if isfile(filename)
    rm(filename)
end

metadata = OrderedDict(
    # Name of the project (SeaDataCloud, SeaDataNet, EMODNET-Chemistry, ...)
    "project" => "SeaDataCloud",

    # URN code for the institution EDMO registry,
    # e.g. SDN:EDMO::1579
    "institution_urn" => "SDN:EDMO::1579",

    # Production group
    "production" => "Diva group. E-mails: a.barth@ulg.ac.be, swatelet@ulg.ac.be",

    # Name and emails from authors
    "Author_e-mail" => ["Your Name1 <name1@example.com>", "Other Name <name2@example.com>"],

    # Source of the observation
    "source" => "observational data from SeaDataNet/EMODNet Chemistry Data Network",

    # Additional comment
    "comment" => "...",

    # SeaDataNet Vocabulary P35 URN
    # http://seadatanet.maris2.nl/v_bodc_vocab_v2/search.asp?lib=p35
    # example: SDN:P35::WATERTEMP
    "parameter_keyword_urn" => "SDN:P35::EPC00001",

    # List of SeaDataNet Parameter Discovery Vocabulary P02 URNs
    # http://seadatanet.maris2.nl/v_bodc_vocab_v2/search.asp?lib=p02
    # example: ["SDN:P02::TEMP"]
    "search_keywords_urn" => ["SDN:P02::PSAL"],

    # List of SeaDataNet Vocabulary C19 area URNs
    # SeaVoX salt and fresh water body gazetteer (C19)
    # http://seadatanet.maris2.nl/v_bodc_vocab_v2/search.asp?lib=C19
    # example: ["SDN:C19::3_1"]
    "area_keywords_urn" => ["SDN:C19::3_3"],

    "product_version" => "1.0",

    # NetCDF CF standard name
    # http://cfconventions.org/Data/cf-standard-names/current/build/cf-standard-name-table.html
    # example "standard_name" = "sea_water_temperature",
    "netcdf_standard_name" => "sea_water_salinity",

    "netcdf_long_name" => "sea water salinity",

    "netcdf_units" => "1e-3",

    # Abstract for the product
    "abstract" => "...",

    # This option provides a place to acknowledge various types of support for the
    # project that produced the data
    "acknowledgment" => "...",

    # Digital Object Identifier of the data product
    "doi" => "...")


@static if VERSION >= v"0.7.0"
    @test_logs (:info,r".*") match_mode=:any DIVAnd.diva3d(
        (lonr,latr,TS),
        (obslon[sel],obslat[sel],obstime[sel]),
        obsvalue[sel],
        (lenx,leny),
        epsilon2,
        filename,varname,
        bathname = bathname,
        bathisglobal = bathisglobal,
        background_len = (lenx,leny),
    )
else
    @test_warn r".*" DIVAnd.diva3d(
        (lonr,latr,TS),
        (obslon[sel],obslat[sel],obstime[sel]),
        obsvalue[sel],
        (lenx,leny),
        epsilon2,
        filename,varname,
        bathname = bathname,
        bathisglobal = bathisglobal,
        background_len = (lenx,leny),
       )
end
