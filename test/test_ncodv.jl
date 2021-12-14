using NCDatasets
using Missings
using DIVAnd
using Test
using DataStructures
using Dates

#=
fnames = [expanduser("~/Downloads/data_from_SDN_2017-11_TS_profiles_non-restricted_med.nc")]
fname = fnames[1]

T = Float64
qv_flags = ["good_value","probably_good_value"];
long_name = "ITS-90 water temperature"

obsvalue,obslon,obslat,obsdepth,obstime,obsids = load(T,fname,long_name; qv_flags = qv_flags)

nothing;


=#

# profiles

fname = tempname()
ds = Dataset(fname, "c")
# Dimensions

ds.dim["N_STATIONS"] = 3
ds.dim["N_SAMPLES"] = 2
ds.dim["STRING2"] = 2

# Declare variables

# Declare variables

nclongitude = defVar(
    ds,
    "longitude",
    Float32,
    ("N_STATIONS",),
    attrib = OrderedDict(
        "long_name" => "Longitude",
        "standard_name" => "longitude",
        "units" => "degrees_east",
        "comment" => "",
        "C_format" => "%.3f",
        "FORTRAN_format" => "F12.3",
        "_FillValue" => Float32(-1.0e10),
    ),
)

nclatitude = defVar(
    ds,
    "latitude",
    Float32,
    ("N_STATIONS",),
    attrib = OrderedDict(
        "long_name" => "Latitude",
        "standard_name" => "latitude",
        "units" => "degrees_north",
        "C_format" => "%.3f",
        "FORTRAN_format" => "F12.3",
        "_FillValue" => Float32(-1.0e10),
    ),
)

ncmetavar4 = defVar(
    ds,
    "metavar4",
    Char,
    ("STRING2", "N_STATIONS"),
    attrib = OrderedDict("long_name" => "LOCAL_CDI_ID"),
)

ncmetavar5 = defVar(
    ds,
    "metavar5",
    Int32,
    ("N_STATIONS",),
    attrib = OrderedDict(
        "long_name" => "EDMO_CODE",
        "C_format" => "%.0f",
        "FORTRAN_format" => "F12.0",
        "_FillValue" => Int32(-2147483646),
    ),
)

ncdate_time = defVar(
    ds,
    "date_time",
    Float64,
    ("N_STATIONS",),
    attrib = OrderedDict(
        "long_name" => "Decimal Gregorian Days of the station",
        "standard_name" => "time",
        "units" => "days since 0190-01-01 00:00:00 UTC",
        "comment" => "Relative Gregorian Days with decimal part",
        "C_format" => "%.5f",
        "FORTRAN_format" => "F12.5",
        "_FillValue" => -1.0e10,
    ),
)

ncvar1 = defVar(
    ds,
    "var1",
    Float32,
    ("N_SAMPLES", "N_STATIONS"),
    attrib = OrderedDict(
        "positive" => "down",
        "long_name" => "Depth",
        "units" => "m",
        "comment" => "Codes: SDN:P01::ADEPZZ01 SDN:P06::ULAA",
        "ancillary_variables" => "var1_qc var1_err",
        "C_format" => "%.2f",
        "FORTRAN_format" => "F12.2",
        "_FillValue" => Float32(-1.0e10),
    ),
)

ncvar1_qc = defVar(
    ds,
    "var1_qc",
    Int8,
    ("N_SAMPLES", "N_STATIONS"),
    attrib = OrderedDict(
        "long_name" => "Quality flag of Depth",
        "standard_name" => "status_flag",
        "comment" => "SEADATANET - SeaDataNet quality codes",
        "flag_values" => Int8[48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 65, 81],
        "flag_meanings" =>
            "no_quality_control good_value probably_good_value probably_bad_value bad_value changed_value value_below_detection value_in_excess interpolated_value missing_value value_phenomenon_uncertain value_below_limit_of_quantification",
        "_FillValue" => Int8(57),
    ),
)

ncvar2 = defVar(
    ds,
    "var2",
    Float32,
    ("N_SAMPLES", "N_STATIONS"),
    attrib = OrderedDict(
        "long_name" => "ITS-90 water temperature",
        "units" => "degrees C",
        "comment" => "Codes: SDN:P35::WATERTEMP SDN:P06::UPAA",
        "ancillary_variables" => "var2_qc var2_err",
        "C_format" => "%.2f",
        "FORTRAN_format" => "F12.2",
        "_FillValue" => Float32(-1.0e10),
    ),
)

ncvar2_qc = defVar(
    ds,
    "var2_qc",
    Int8,
    ("N_SAMPLES", "N_STATIONS"),
    attrib = OrderedDict(
        "long_name" => "Quality flag of ITS-90 water temperature",
        "standard_name" => "status_flag",
        "comment" => "SEADATANET - SeaDataNet quality codes",
        "flag_values" => Int8[48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 65, 81],
        "flag_meanings" =>
            "no_quality_control good_value probably_good_value probably_bad_value bad_value changed_value value_below_detection value_in_excess interpolated_value missing_value value_phenomenon_uncertain value_below_limit_of_quantification",
        "_FillValue" => Int8(57),
    ),
)


# Global attributes

ds.attrib["Conventions"] = "CF-1.7"
ds.attrib["comment"] = "ODV NetCDF Export File V2.0"
ds.attrib["Software"] = "Ocean Data View 5.1.4 - 64 bit (Linux)"
ds.attrib["DataField"] = "Ocean"
ds.attrib["DataType"] = "Profiles"

# Define variables

nclongitude[:] = [1, 2, 3]
nclatitude[:] = [1, 2, 3]
ncmetavar4[:] = ['a' 'b' 'c'; 'x' 'y' 'z']
ncmetavar5[:] = [111, 222, 333]
ncdate_time[:] = [DateTime(2000, 1, 1), DateTime(2000, 1, 2), DateTime(2000, 1, 3)]
ncvar1[:] = [0.0 0.0 0.0; 1.0 1.0 1.0]
ncvar1_qc[:] = [49 49 49; 49 49 48]
ncvar2[:] = [10.0 11.0 12.0; 20.0 21.0 22.0]
ncvar2_qc[:] = [49 48 49; 49 49 49]
close(ds)

T = Float64
qv_flags = ["good_value", "probably_good_value"];
long_name = "ITS-90 water temperature"

obsvalue, obslon, obslat, obsdepth, obstime, obsids =
    NCODV.load(T, fname, long_name; qv_flags = qv_flags)

@test obsvalue == [10.0, 20.0, 21.0, 12.0]
@test obslon == [1.0, 1.0, 2.0, 3.0]
@test obslat == [1.0, 1.0, 2.0, 3.0]
@test obstime == [
    DateTime(2000, 1, 1),
    DateTime(2000, 1, 1),
    DateTime(2000, 1, 2),
    DateTime(2000, 1, 3),
]
@test obsdepth == [0.0, 1.0, 1.0, 0.0]
@test obsids == ["111-ax", "111-ax", "222-by", "333-cz"]


# time series

using NCDatasets, DataStructures
fname_TS = tempname()


ds = NCDataset(fname_TS,"c", attrib = OrderedDict(
    "Conventions"               => "CF-1.7",
    "comment"                   => "ODV NetCDF Export File V2.0",
    "CreateTime"                => "2021-03-15T21:32:19",
    "Software"                  => "Ocean Data View 5.3.11 - 64 bit (Windows)",
    "SourceLastModified"        => "2021-03-15T15:45:55",
    "DataField"                 => "Ocean",
    "DataType"                  => "Profiles",
    "Description"               => "",
    "featureType"               => "profile",
))

# Dimensions

ds.dim["N_STATIONS"] = 3
ds.dim["N_SAMPLES"] = 2
ds.dim["STRING2"] = 2


# Declare variables


nclongitude = defVar(ds,"longitude", Float32, ("N_STATIONS",), attrib = OrderedDict(
    "long_name"                 => "Longitude",
    "standard_name"             => "longitude",
    "units"                     => "degrees_east",
    "comment"                   => "",
    "C_format"                  => "%.3f",
    "FORTRAN_format"            => "F12.3",
    "_FillValue"                => Float32(-1.0e10),
))

nclatitude = defVar(ds,"latitude", Float32, ("N_STATIONS",), attrib = OrderedDict(
    "long_name"                 => "Latitude",
    "standard_name"             => "latitude",
    "units"                     => "degrees_north",
    "comment"                   => "",
    "C_format"                  => "%.3f",
    "FORTRAN_format"            => "F12.3",
    "_FillValue"                => Float32(-1.0e10),
))

ncmetavar4 = defVar(ds,"metavar4", Char, ("STRING2", "N_STATIONS"), attrib = OrderedDict(
    "long_name"                 => "LOCAL_CDI_ID",
    "units"                     => "",
    "comment"                   => "",
))

ncmetavar5 = defVar(ds,"metavar5", Int32, ("N_STATIONS",), attrib = OrderedDict(
    "long_name"                 => "EDMO_code",
    "units"                     => "",
    "comment"                   => "",
    "C_format"                  => "%.0f",
    "FORTRAN_format"            => "F12.0",
    "_FillValue"                => Int32(-2147483646),
))


ncdate_time = defVar(ds,"date_time", Float64, ("N_STATIONS",), attrib = OrderedDict(
    "long_name"                 => "Decimal Gregorian Days of the station",
    "standard_name"             => "time",
    "units"                     => "days since 1974-01-01 00:00:00 UTC",
    "comment"                   => "Relative Gregorian Days with decimal part",
    "C_format"                  => "%.5f",
    "FORTRAN_format"            => "F12.5",
    "_FillValue"                => -1.0e10,
))

ncvar1 = defVar(ds,"var1", Float32, ("N_SAMPLES", "N_STATIONS"), attrib = OrderedDict(
    "long_name"                 => "time_ISO8601",
    "standard_name"             => "time",
    "units"                     => "years since 0000-01-01",
    "comment"                   => "Codes: SDN:P01::DTUT8601 SDN:P06::TISO",
    "ancillary_variables"       => "var1_qc",
    "C_format"                  => "%.2f",
    "FORTRAN_format"            => "F12.2",
    "_FillValue"                => Float32(-1.0e10),
))

ncvar1_qc = defVar(ds,"var1_qc", Int8, ("N_SAMPLES", "N_STATIONS"), attrib = OrderedDict(
    "long_name"                 => "Quality flag of time_ISO8601",
    "standard_name"             => "status_flag",
    "comment"                   => "SEADATANET - SeaDataNet quality codes",
    "flag_values"               => Int8[48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 65, 81],
    "flag_meanings"             => "no_quality_control good_value probably_good_value probably_bad_value bad_value changed_value value_below_detection value_in_excess interpolated_value missing_value value_phenomenon_uncertain value_below_limit_of_quantification",
    "_FillValue"                => Int8(57),
))

ncvar2 = defVar(ds,"var2", Float32, ("N_SAMPLES", "N_STATIONS"), attrib = OrderedDict(
    "long_name"                 => "Depth",
    "units"                     => "m",
    "comment"                   => "Codes: SDN:P01::ADEPZZ01 SDN:P06::ULAA | SDN:P01::DEPHPR01 SDN:P06::ULAA | SDN:P01::DEPHWO01 SDN:P06::ULAA",
    "ancillary_variables"       => "var2_qc",
    "C_format"                  => "%.2f",
    "FORTRAN_format"            => "F12.2",
    "_FillValue"                => Float32(-1.0e10),
))

ncvar2_qc = defVar(ds,"var2_qc", Int8, ("N_SAMPLES", "N_STATIONS"), attrib = OrderedDict(
    "long_name"                 => "Quality flag of Depth",
    "standard_name"             => "status_flag",
    "comment"                   => "SEADATANET - SeaDataNet quality codes",
    "flag_values"               => Int8[48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 65, 81],
    "flag_meanings"             => "no_quality_control good_value probably_good_value probably_bad_value bad_value changed_value value_below_detection value_in_excess interpolated_value missing_value value_phenomenon_uncertain value_below_limit_of_quantification",
    "_FillValue"                => Int8(57),
))


ncvar3 = defVar(ds,"var3", Float32, ("N_SAMPLES", "N_STATIONS"), attrib = OrderedDict(
    "long_name"                 => "Water body dissolved oxygen concentration",
    "units"                     => "umol/l",
    "comment"                   => "Codes: SDN:P35::EPC00002 SDN:P06::UPOX",
    "ancillary_variables"       => "var3_qc",
    "C_format"                  => "%.2f",
    "FORTRAN_format"            => "F12.2",
    "_FillValue"                => Float32(-1.0e10),
))

ncvar3_qc = defVar(ds,"var3_qc", Int8, ("N_SAMPLES", "N_STATIONS"), attrib = OrderedDict(
    "long_name"                 => "Quality flag of Water body dissolved oxygen concentration",
    "standard_name"             => "status_flag",
    "comment"                   => "SEADATANET - SeaDataNet quality codes",
    "flag_values"               => Int8[48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 65, 81],
    "flag_meanings"             => "no_quality_control good_value probably_good_value probably_bad_value bad_value changed_value value_below_detection value_in_excess interpolated_value missing_value value_phenomenon_uncertain value_below_limit_of_quantification",
    "_FillValue"                => Int8(57),
))

# Define variables

nclongitude[:] = [1, 2, 3]
nclatitude[:] = [1, 2, 3]
ncmetavar4[:] = ['a' 'b' 'c'; 'x' 'y' 'z']
ncmetavar5[:] = [111, 222, 333]
ncdate_time[:] = [DateTime(1900, 1, 1), DateTime(1900, 1, 2), DateTime(1900, 1, 3)]

# time
ncvar1[:] =    [2000.0 2000.0 2000.0; 2001.0 2001.0 2001.0]
ncvar1_qc[:] = [  49     49     49;     49     49     48  ]

# depth
ncvar2[:] =    [  10.0   11.0   12.0;  20.0    21.0   22.0]
ncvar2_qc[:] = [  49     48     49;    49      49     49  ]

# measurement
ncvar3[:] =    [  10.0   11.0   12.0;  20.0    21.0   22.0]
ncvar3_qc[:] = [  49     49     49;    49      49     48  ]

close(ds)


# load file
T = Float64
long_name = "Water body dissolved oxygen concentration"

obsvalue,obslon,obslat,obsdepth,obstime,obsids = NCODV.load(T,fname_TS,long_name)

almostgood(qf) = !ismissing(qf) && ((qf == 49) || (qf == 50))

nc = NCDataset(fname_TS)
ox = nc["var3"][:];
ox_qf = nc["var3_qc"][:];

z = nc["var2"][:];
z_qf = nc["var2_qc"][:];

time = nc["var1_qc"].var[:];
time_qf = nc["var1_qc"][:];

good = almostgood.(ox_qf) .& .!ismissing.(ox) .& almostgood.(z_qf) .& .!ismissing.(z) .& almostgood.(time_qf) .& .!ismissing.(time)

@test ox[good] == obsvalue


# Test warning of missing flags
fname = "ncodv_without_flags.nc"
if !isfile(fname)
    download("https://dox.ulg.ac.be/index.php/s/nGmsH5ydCAjgft5/download",fname)
end
T = Float32
long_name = "Water body salinity"

@test_logs (:warn, r".*status_flags.*") match_mode = :any NCODV.load(
   T,fname,long_name)
