using NCDatasets
using Missings
using DIVAnd
if VERSION >= v"0.7.0-beta.0"
    using Test
else
    using Base.Test
end





#=
fnames = [expanduser("~/Downloads/data_from_SDN_2017-11_TS_profiles_non-restricted_med.nc")]
fname = fnames[1]

T = Float64
qv_flags = ["good_value","probably_good_value"];
long_name = "ITS-90 water temperature"

obsvalue,obslon,obslat,obsdepth,obstime,obsids = load(T,fname,long_name; qv_flags = qv_flags)

nothing;


=#

fname = tempname()
ds = Dataset(fname,"c")
# Dimensions

ds.dim["N_STATIONS"] = 3
ds.dim["N_SAMPLES"] = 2
ds.dim["STRING2"] = 2

# Declare variables


nclongitude = defVar(ds,"longitude", Float32, ("N_STATIONS",))
nclongitude.attrib["long_name"] = "Longitude"
nclongitude.attrib["standard_name"] = "longitude"
nclongitude.attrib["units"] = "degrees_east"
nclongitude.attrib["comment"] = ""
nclongitude.attrib["C_format"] = "%.3f"
nclongitude.attrib["FORTRAN_format"] = "F12.3"
nclongitude.attrib["_FillValue"] = Float32(-1.0e10)

nclatitude = defVar(ds,"latitude", Float32, ("N_STATIONS",))
nclatitude.attrib["long_name"] = "Latitude"
nclatitude.attrib["standard_name"] = "latitude"
nclatitude.attrib["units"] = "degrees_north"
nclatitude.attrib["C_format"] = "%.3f"
nclatitude.attrib["FORTRAN_format"] = "F12.3"
nclatitude.attrib["_FillValue"] = Float32(-1.0e10)

ncmetavar4 = defVar(ds,"metavar4", Char, ("STRING2", "N_STATIONS"))
ncmetavar4.attrib["long_name"] = "LOCAL_CDI_ID"

ncmetavar5 = defVar(ds,"metavar5", Int32, ("N_STATIONS",))
ncmetavar5.attrib["long_name"] = "EDMO_CODE"
ncmetavar5.attrib["C_format"] = "%.0f"
ncmetavar5.attrib["FORTRAN_format"] = "F12.0"
ncmetavar5.attrib["_FillValue"] = -2147483646

ncdate_time = defVar(ds,"date_time", Float64, ("N_STATIONS",))
ncdate_time.attrib["long_name"] = "Decimal Gregorian Days of the station"
ncdate_time.attrib["standard_name"] = "time"
ncdate_time.attrib["units"] = "days since 0190-01-01 00:00:00 UTC"
ncdate_time.attrib["comment"] = "Relative Gregorian Days with decimal part"
ncdate_time.attrib["C_format"] = "%.5f"
ncdate_time.attrib["FORTRAN_format"] = "F12.5"
ncdate_time.attrib["_FillValue"] = -1.0e10

ncvar1 = defVar(ds,"var1", Float32, ("N_SAMPLES", "N_STATIONS"))
ncvar1.attrib["positive"] = "down"
ncvar1.attrib["long_name"] = "Depth"
ncvar1.attrib["units"] = "m"
ncvar1.attrib["comment"] = "Codes: SDN:P01::ADEPZZ01 SDN:P06::ULAA"
ncvar1.attrib["ancillary_variables"] = "var1_qc var1_err"
ncvar1.attrib["C_format"] = "%.2f"
ncvar1.attrib["FORTRAN_format"] = "F12.2"
ncvar1.attrib["_FillValue"] = Float32(-1.0e10)

ncvar1_qc = defVar(ds,"var1_qc", Int8, ("N_SAMPLES", "N_STATIONS"))
ncvar1_qc.attrib["long_name"] = "Quality flag of Depth"
ncvar1_qc.attrib["standard_name"] = "status_flag"
ncvar1_qc.attrib["comment"] = "SEADATANET - SeaDataNet quality codes"
ncvar1_qc.attrib["flag_values"] = Int8[48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 65, 81]
ncvar1_qc.attrib["flag_meanings"] = "no_quality_control good_value probably_good_value probably_bad_value bad_value changed_value value_below_detection value_in_excess interpolated_value missing_value value_phenomenon_uncertain value_below_limit_of_quantification"
ncvar1_qc.attrib["_FillValue"] = Int8(57)


ncvar2 = defVar(ds,"var2", Float32, ("N_SAMPLES", "N_STATIONS"))
ncvar2.attrib["long_name"] = "ITS-90 water temperature"
ncvar2.attrib["units"] = "degrees C"
ncvar2.attrib["comment"] = "Codes: SDN:P35::WATERTEMP SDN:P06::UPAA"
ncvar2.attrib["ancillary_variables"] = "var2_qc var2_err"
ncvar2.attrib["C_format"] = "%.2f"
ncvar2.attrib["FORTRAN_format"] = "F12.2"
ncvar2.attrib["_FillValue"] = Float32(-1.0e10)

ncvar2_qc = defVar(ds,"var2_qc", Int8, ("N_SAMPLES", "N_STATIONS"))
ncvar2_qc.attrib["long_name"] = "Quality flag of ITS-90 water temperature"
ncvar2_qc.attrib["standard_name"] = "status_flag"
ncvar2_qc.attrib["comment"] = "SEADATANET - SeaDataNet quality codes"
ncvar2_qc.attrib["flag_values"] = Int8[48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 65, 81]
ncvar2_qc.attrib["flag_meanings"] = "no_quality_control good_value probably_good_value probably_bad_value bad_value changed_value value_below_detection value_in_excess interpolated_value missing_value value_phenomenon_uncertain value_below_limit_of_quantification"
ncvar2_qc.attrib["_FillValue"] = Int8(57)

# Global attributes

ds.attrib["Conventions"] = "CF-1.7"
ds.attrib["comment"] = "ODV NetCDF Export File V2.0"
ds.attrib["Software"] = "Ocean Data View 5.1.4 - 64 bit (Linux)"
ds.attrib["DataField"] = "Ocean"
ds.attrib["DataType"] = "Profiles"

# Define variables

nclongitude[:] = [1,2,3]
nclatitude[:] = [1,2,3]
ncmetavar4[:] = ['a' 'b' 'c'; 'x' 'y' 'z']
ncmetavar5[:] = [111,222,333]
ncdate_time[:] = [DateTime(2000,1,1),DateTime(2000,1,2),DateTime(2000,1,3)]
ncvar1[:] = [0. 0. 0.; 1. 1. 1.]
ncvar1_qc[:] = [49 49 49; 49 49 48]
ncvar2[:] =  [10. 11. 12.; 20. 21. 22.]
ncvar2_qc[:] = [49 48 49; 49 49 49]
close(ds)

T = Float64
qv_flags = ["good_value","probably_good_value"];
long_name = "ITS-90 water temperature"

obsvalue,obslon,obslat,obsdepth,obstime,obsids = NCODV.load(T,fname,long_name; qv_flags = qv_flags)

@test obsvalue == [10.0, 20.0, 21.0, 12.0]
@test obslon == [1.0, 1.0, 2.0, 3.0]
@test obslat == [1.0, 1.0, 2.0, 3.0]
@test obstime == [DateTime(2000,1,1),DateTime(2000,1,1),DateTime(2000,1,2),DateTime(2000,1,3)]
@test obsdepth == [0.0, 1.0, 1.0, 0.0]
@test obsids == ["111-ax", "111-ax", "222-by", "333-cz"]
