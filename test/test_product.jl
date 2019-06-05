import DIVAnd
if VERSION >= v"0.7.0-beta.0"
    using Test
    using DelimitedFiles
else
    using Base.Test
    using Compat
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

obsvalue,obslon,obslat,obsdepth,obstime,obsids = DIVAnd.loadobs(Float64,obsname,"Salinity")

# for testing only
obsids[1] = "100-123"
obsids[2] = "100-124"
obsids[3] = "999-missing"
obsids[4:end] .= "101-125"


dx = dy = 0.5
lonr = 3:dx:11.8
latr = 42.0:dy:44.0
depthr = [0.,20.,30.]
epsilon2 = 0.01

# put one point on land
index_land_point = 897297
obslat[index_land_point] = 43.6333
obslon[index_land_point] = 6

# put one to NaN
index_NaN = 897298
obsvalue[index_NaN] = NaN

# put an outlier
index_outlier = 897299
obsvalue[index_outlier] = 50.


surfextend = true
sz = (length(lonr),length(latr),length(depthr))

lenx = fill(200_000,sz)
leny = fill(200_000,sz)
lenz = [10+depthr[k]/15 for i = 1:sz[1], j = 1:sz[2], k = 1:sz[3]]

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
filename = "Water_body_$(replace(varname,' ' => '_')).4Danl.nc"


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


# edit the bathymetry
mask,(pm,pn,po),(xi,yi,zi) = DIVAnd.domain(bathname,bathisglobal,lonr,latr,depthr)
mask[3,3,1] = false

ncglobalattrib,ncvarattrib = DIVAnd.SDNMetadata(metadata,filename,varname,lonr,latr)

if isfile(filename)
   rm(filename) # delete the previous analysis
end

dbinfo = @static if VERSION >= v"0.7.0"
    @test_logs (:info,r".*netCDF*") match_mode=:any DIVAnd.diva3d(
        (lonr,latr,depthr,TS),
        (obslon,obslat,obsdepth,obstime),
        obsvalue,
        (lenx,leny,lenz),
        epsilon2,
        filename,varname,
        bathname = bathname,
        bathisglobal = bathisglobal,
        ncvarattrib = ncvarattrib,
        ncglobalattrib = ncglobalattrib,
        mask = mask,
        surfextend = surfextend,
    )
else
    @test_warn r".netCDF.*" DIVAnd.diva3d(
        (lonr,latr,depthr,TS),
        (obslon,obslat,obsdepth,obstime),
        obsvalue,
        (lenx,leny,lenz),
        epsilon2,
        filename,varname,
        bathname = bathname,
        bathisglobal = bathisglobal,
        ncvarattrib = ncvarattrib,
        ncglobalattrib = ncglobalattrib,
        mask = mask,
        surfextend = surfextend,
    )
end

obsused = dbinfo[:used]
#DIVAnd.saveobs(filename,(obslon,obslat,obsdepth,obstime),obsids,used = obsused)

DIVAnd.saveobs(filename,(obslon,obslat,obsdepth,obstime),obsids)


project = "SeaDataCloud"
xmlfilename = "test.xml"
ignore_errors = true

additionalcontacts = [
    DIVAnd.getedmoinfo(1977,"originator"),
    DIVAnd.getedmoinfo(4630,"originator"),
]

@static if VERSION >= v"0.7.0"
    @test_logs (:info,r".*") match_mode=:any DIVAnd.divadoxml(
        filename,varname,project,cdilist,xmlfilename,
        ignore_errors = ignore_errors,
        additionalcontacts = additionalcontacts
    )
else
    @test_warn r".*" DIVAnd.divadoxml(
        filename,varname,project,cdilist,xmlfilename,
        ignore_errors = ignore_errors,
        additionalcontacts = additionalcontacts
    )
end

errname = "$(replace(filename,r"\.nc$" => "")).cdi_import_errors_test.csv"

errdata,header = readdlm(errname,'\t'; header = true)

# check if the missing CDI was identified
@test errdata[1] == 999
@test errdata[2] == "missing"

# check if editing of the mask was successful
@test ismissing(Dataset(filename)["Salinity"][3,3,1,1])

xmlstr = read(xmlfilename,String);

keyword_code = split(metadata["parameter_keyword_urn"],':')[end]
@test occursin(keyword_code,xmlstr)
@test occursin("1977",xmlstr)
@test occursin("4630",xmlstr)

# new analysis with background from file

filename2 = "Water_body_$(replace(varname,' ' => '_'))2.4Danl.nc"
if isfile(filename2)
   rm(filename2) # delete the previous analysis
end


dbinfo =
    @static if VERSION >= v"0.7.0"
        @test_logs (:info,r".*") (:warn,r".*Be patient.*") match_mode=:any DIVAnd.diva3d(
            (lonr,latr,depthr,TS),
            (obslon,obslat,obsdepth,obstime),
            obsvalue,
            (),
            epsilon2,
            filename2,varname,
            bathname = bathname,
            bathisglobal = bathisglobal,
            ncvarattrib = ncvarattrib,
            ncglobalattrib = ncglobalattrib,
            background = DIVAnd.backgroundfile(filename,varname),
            fitcorrlen = true,
            background_len = (lenx,leny,lenz),
            fithorz_param = Dict(
                :maxnsamp => 500,
                :epsilon2 => ones(size(obsvalue))
            ),
            fitvert_param = Dict(
                :maxnsamp => 100,
            ),
            mask = mask,
            niter_e = 2,
            QCMETHOD = 0,
            surfextend = surfextend,
        )
    else
        @test_warn r".*Be patient.*" DIVAnd.diva3d(
            (lonr,latr,depthr,TS),
            (obslon,obslat,obsdepth,obstime),
            obsvalue,
            (),
            epsilon2,
            filename2,varname,
            bathname = bathname,
            bathisglobal = bathisglobal,
            ncvarattrib = ncvarattrib,
            ncglobalattrib = ncglobalattrib,
            background = DIVAnd.backgroundfile(filename,varname),
            fitcorrlen = true,
            background_len = (lenx,leny,lenz),
            fithorz_param = Dict(
                :maxnsamp => 500,
            ),
            fitvert_param = Dict(
                :maxnsamp => 100,
            ),
            mask = mask,
            niter_e = 2,
            QCMETHOD = 0,
            surfextend = surfextend,
        )
    end

qcvalue = dbinfo[:qcvalues]
used = dbinfo[:used]
residuals = dbinfo[:residuals]

@test isnan(residuals[index_land_point])
@test isnan(residuals[index_NaN])
@test !used[index_NaN]
@test all(isfinite.(qcvalue[used]))
@test qcvalue[index_outlier] > 9

nothing

