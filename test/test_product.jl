import DIVAnd
using DelimitedFiles
using DataStructures
using Interpolations
using Missings
using NCDatasets
using Test
using Interpolations

varname = "Salinity"
filename = "WOD-Salinity.nc"


bathname = joinpath(
    dirname(@__FILE__),
    "..",
    "..",
    "DIVAnd-example-data",
    "Global",
    "Bathymetry",
    "gebco_30sec_16.nc",
)
bathisglobal = true

obsname = joinpath(
    dirname(@__FILE__),
    "..",
    "..",
    "DIVAnd-example-data",
    "Provencal",
    "WOD-Salinity.nc",
)

cdilist = joinpath(dirname(@__FILE__), "..", "data", "CDI-list-export.csv")


if !isfile(bathname)
    @info("download bathymetry $bathname")
    bathname = download("https://dox.ulg.ac.be/index.php/s/U0pqyXhcQrXjEUX/download")
end


if !isfile(obsname)
    @info("download observations $obsname")
    obsname = download("https://dox.ulg.ac.be/index.php/s/PztJfSEnc8Cr3XN/download")
end

obsvalue, obslon, obslat, obsdepth, obstime, obsids =
    DIVAnd.loadobs(Float64, obsname, "Salinity")

# for testing only
obsids[1] = "100-123"
obsids[2] = "100-124"
obsids[3] = "999-missing"
obsids[4:end] .= "101-125"


dx = dy = 0.5
lonr = 3:dx:11.8
latr = 42.0:dy:44.0
depthr = [0.0, 20.0, 30.0]
depthr = [0.0, 20.0, 30.0, 500, 1000, 2000]
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
obsvalue[index_outlier] = 50.0


surfextend = true
sz = (length(lonr), length(latr), length(depthr))

lenx = fill(200_000, sz)
leny = fill(200_000, sz)
lenz = [10 + depthr[k] / 15 for i = 1:sz[1], j = 1:sz[2], k = 1:sz[3]]

years = 1993:1993

year_window = 10

# winter: January-March    1,2,3
# spring: April-June       4,5,6
# summer: July-September   7,8,9
# autumn: October-December 10,11,12

monthlists = [[1, 2, 3], [4, 5, 6], [7, 8, 9], [10, 11, 12]];


TS = DIVAnd.TimeSelectorYW(years, year_window, monthlists)

varname = "Salinity"

# File name
filename = tempname()

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
    "doi" => "...",
)


# edit the bathymetry
mask, (pm, pn, po), (xi, yi, zi) = DIVAnd.domain(bathname, bathisglobal, lonr, latr, depthr)
mask[3, 3, 1] = false

ncglobalattrib, ncvarattrib = DIVAnd.SDNMetadata(metadata, filename, varname, lonr, latr)

error_thresholds = [("L1", 0.3), ("L2", 0.5), ("L3", 0.1)]

if isfile(filename)
    rm(filename) # delete the previous analysis
end

dbinfo = @test_logs (:info, r".*netCDF*") match_mode = :any DIVAnd.diva3d(
    (lonr, latr, depthr, TS),
    (obslon, obslat, obsdepth, obstime),
    obsvalue,
    (lenx, leny, lenz),
    epsilon2,
    filename,
    varname,
    bathname = bathname,
    bathisglobal = bathisglobal,
    ncvarattrib = ncvarattrib,
    ncglobalattrib = ncglobalattrib,
    mask = mask,
    surfextend = surfextend,
    stat_per_timeslice = true,
    error_thresholds = error_thresholds,
)



# observation statistics in data bins
meanv,count = DIVAnd.binning((lonr,latr,depthr),(obslon,obslat,obsdepth),obsvalue)

@test sum(count) <= length(obsvalue)

meanv,count = DIVAnd.binning((lonr,latr,depthr,TS),(obslon,obslat,obsdepth,obstime),obsvalue)

@test sum(count) <= length(obsvalue)
@test maximum(filter(isfinite,meanv)) <= maximum(filter(isfinite,obsvalue))
@test minimum(filter(isfinite,meanv)) >= minimum(filter(isfinite,obsvalue))

isoutlier = falses(size(obsvalue))

xyi = (lonr,latr,depthr,TS)
xy = (obslon, obslat, obsdepth, obstime)
DIVAnd.saveobsstat(filename,xyi,xy; isoutlier = isoutlier)


# save observations
obsused = dbinfo[:used]
#DIVAnd.saveobs(filename,(obslon,obslat,obsdepth,obstime),obsids,used = obsused)
DIVAnd.saveobs(filename, (obslon, obslat, obsdepth, obstime), obsids)

# derived parameters
filename2 = tempname()
DIVAnd.derived(filename,varname,filename2,error_thresholds = error_thresholds)

# cutting results
filename_cut = tempname()
polygon_lon = lonr[[5, end, end, 5]]
polygon_lat = latr[[1, 1, end, end]]

maskkeep = DIVAnd.inpolygon(polygon_lon,polygon_lat,lonr,latr)

DIVAnd.cut(filename2,varname,filename_cut,polygon_lon,polygon_lat)



project = "SeaDataCloud"
xmlfilename = tempname()
ignore_errors = true

additionalcontacts =
    [DIVAnd.getedmoinfo(1977, "originator"), DIVAnd.getedmoinfo(4630, "originator")]

errname = split(filename, ".nc")[1] * ".cdi_import_errors.csv"

@test_logs (:info, r".*") match_mode = :any DIVAnd.divadoxml(
    filename,
    varname,
    project,
    cdilist,
    xmlfilename,
    ignore_errors = ignore_errors,
    additionalcontacts = additionalcontacts,
)

errdata, header = readdlm(errname, '\t'; header = true)

# check if the missing CDI was identified
@test errdata[1] == 999
@test errdata[2] == "missing"

# check if editing of the mask was successful
ds = Dataset(filename)
@test ismissing(ds["Salinity"][3, 3, 1, 1])
@test haskey(ds,"Salinity_L3")

close(ds)

xmlstr = read(xmlfilename, String);

keyword_code = split(metadata["parameter_keyword_urn"], ':')[end]
@test occursin(keyword_code, xmlstr)
@test occursin("1977", xmlstr)
@test occursin("4630", xmlstr)

# new analysis with background from file

filename2 = tempname()
if isfile(filename2)
    rm(filename2) # delete the previous analysis
end


dbinfo = @test_logs (:info, r".*") match_mode = :any DIVAnd.diva3d(
    (lonr, latr, depthr, TS),
    (obslon, obslat, obsdepth, obstime),
    obsvalue,
    (),
    epsilon2,
    filename2,
    varname,
    bathname = bathname,
    bathisglobal = bathisglobal,
    ncvarattrib = ncvarattrib,
    ncglobalattrib = ncglobalattrib,
    background = DIVAnd.backgroundfile(filename, varname),
    fitcorrlen = true,
    background_len = (lenx, leny, lenz),
    fithorz_param = Dict(:maxnsamp => 500, :epsilon2 => ones(size(obsvalue))),
    fitvert_param = Dict(:maxnsamp => 100),
    mask = mask,
    niter_e = 2,
    QCMETHOD = 0,
    surfextend = surfextend,
    stat_per_timeslice = true,
)

qcvalue = dbinfo[:qcvalues]
used = dbinfo[:used]
residuals = dbinfo[:residuals]

@test isnan(residuals[index_land_point])
@test isnan(residuals[index_NaN])
@test !used[index_NaN]
@test all(isfinite.(qcvalue[used]))
@test qcvalue[index_outlier] > 9

rm(xmlfilename)
rm(errname)
rm(filename2)

# ------------------
# interpolate background from a NetCDF file

# reuse previously created file filename
varname = "Salinity"

ds = Dataset(filename)
lon = nomissing(ds["lon"][:])
lat = nomissing(ds["lat"][:])
depth = nomissing(ds["depth"][:])
time = nomissing(ds["time"][:])

v = ds["Salinity"][:, :, :, :]
close(ds)

i = 3
j = 2
k = 2
n = 2

loni = [lon[i]]
lati = [lat[j]]
depthi = [10.0]
timei = [time[n]]


x = (lon, lat, depth)
xi = (loni, lati, depthi)

vn = zeros(size(v[:, :, :, n]))
vn[:] = map((x -> ismissing(x) ? NaN : x), v[:, :, :, n]);

fi = DIVAnd.interp(x, vn, xi)


firef = [(v[i, j, 1, n] + v[i, j, 2, n]) / 2]
@test fi ≈ firef

background = DIVAnd.backgroundfile(filename, varname)
vn2, fi = background(xi, n, firef, DIVAnd.Anam.notransform()[1])

@test fi ≈ [0] atol = 1e-5


#removing the file creates issues on Windows
#rm(filename)

nothing
