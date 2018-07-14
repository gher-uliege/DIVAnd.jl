import DIVAnd
using DataStructures

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
    "parameter_keyword_urn" => "SDN:P35::WATERTEMP",

    # List of SeaDataNet Parameter Discovery Vocabulary P02 URNs
    # http://seadatanet.maris2.nl/v_bodc_vocab_v2/search.asp?lib=p02
    # example: ["SDN:P02::TEMP"]
    "search_keywords_urn" => ["SDN:P02::TEMP"],

    # List of SeaDataNet Vocabulary C19 area URNs
    # SeaVoX salt and fresh water body gazetteer (C19)
    # http://seadatanet.maris2.nl/v_bodc_vocab_v2/search.asp?lib=C19
    # example: ["SDN:C19::3_1"]
    "area_keywords_urn" => ["SDN:C19::3_1"],

    "product_version" => "1.0",

    # NetCDF CF standard name
    # http://cfconventions.org/Data/cf-standard-names/current/build/cf-standard-name-table.html
    # example "standard_name" = "sea_water_temperature",
    "netcdf_standard_name" => "sea_water_temperature",

    "netcdf_long_name" => "sea water temperature",

    "netcdf_units" => "degree Celsius",

    # Abstract for the product
    "abstract" => "...",

    # This option provides a place to acknowledge various types of support for the
    # project that produced the data
    "acknowledgment" => "...",

    # Digital Object Identifier of the data product
    "doi" => "...")

filename = "test.nc"
varname = "temp"
ncglobalattrib,ncvarattrib = DIVAnd.SDNMetadata(metadata,filename,varname,1:10,1:10;
                                         field = ones(10,10))

@test ncvarattrib["units"] == "degree Celsius"
@test typeof(ncglobalattrib["area_keywords"]) == String

ncglobalattrib,ncvarattrib = DIVAnd.SDNMetadata(metadata,filename,varname,1:10,1:10;
                                         default_field_min = -123,
                                         default_field_max = 123)

@test occursin("123",ncglobalattrib["preview"])


