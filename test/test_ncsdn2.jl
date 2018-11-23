using NCDatasets
using Missings
import DIVAnd
if VERSION >= v"0.7.0-beta.0"
    using Test
else
    using Base.Test
end

NCSDN = DIVAnd.NCSDN



function varbyattrib_first(ds; kwargs...)
    vs = varbyattrib(ds; kwargs...)
    if length(vs) == 0
        str = join(["attribute '$k' equal to '$v'" for (k,v) in kwargs]," and ")
        error("No NetCDF variable found with $(str) in $(path(ds))")
    end
    return vs[1]
end


# # files always hava variable with the long_name  "LOCAL_CDI_ID" and "EDMO_CODE" (all upper-case)
# # long_name for the primary variable to analysis are always P35 names
# # longitude, latitude and time (including dates) have the standard attribute "longitude", "latitude" and "time" respectively


@inline function load!(ncvar::NCDatasets.Variable{T,2}, data, i::Colon,j::UnitRange) where T
    n_samples = size(ncvar,1)
    # reversed and 0-based
    NCDatasets.nc_get_vara!(ncvar.ncid,ncvar.varid,[first(j)-1,0],[length(j),n_samples],data)
#    @show datac[1:10],size(datac)
end


@inline function load!(ncvar::NCDatasets.Variable{T,2}, data, indices::Union{Integer, UnitRange, Colon}...) where T
    ind = to_indices(ncvar,indices)
    start = [first(ind[2])-1,first(ind[1])-1]
    count = [length(ind[2]),length(ind[1])]
    stride = [step(ind[2]),step(ind[1])]
    NCDatasets.nc_get_vars!(ncvar.ncid,ncvar.varid,start,count,stride,data)
end

function loadmis(ncvar::NCDatasets.Variable{T,2},fillval) where T
    n_samples,n_stations = size(ncvar)
    n_stations = 100000
    data = Vector{Vector{T}}(undef,n_stations);
    tmp = Array{T,2}(undef,(n_samples,1))

    @inbounds for i = 1:n_stations
        # reversed and 0-based
        NCDatasets.nc_get_vars!(ncvar.ncid,ncvar.varid,[i-1,0],[1,n_samples],[1,1],tmp)
        #data[i] = tmp[tmp .!= fillval]
        data[i] = filter(x -> x != fillval,tmp)
    end
    return data
end

#=
reference for n_stations = 10000
1.467 s (14975 allocations: 3.78 MiB)
1.459 s vara
1.463 s (19975 allocations: 4.03 MiB) - load!
=#
function loadmis2(ncvar::NCDatasets.Variable{T,2},fillval) where T
    nchunk = 10
    n_samples,n_stations = size(ncvar)
#    n_stations = 100000
    n_stations = 10000
    data = Vector{Vector{T}}(undef,n_stations);
    data_chunk = Array{T,2}(undef,(n_samples,nchunk))

    profile = Vector{T}(undef,n_samples)

    @inbounds for i = 1:nchunk:n_stations
        if (i-1) % (nchunk*1000) == 0
            println("$(i-1) out of $n_stations - $(100*(i-1)/n_stations) %")
        end
        nc = i:min(i+nchunk-1,n_stations)
        clen = length(nc)

        #NCDatasets.nc_get_vara(ncvar.ncid,ncvar.varid,[i-1,0],[clen,n_samples],data_chunk)
        load!(ncvar,data_chunk,:,nc)

        for k = 1:clen
            iprofile = 0
            for l = 1:n_samples
                if data_chunk[l,k] != fillval
                    iprofile = iprofile+1
                    profile[iprofile] = data_chunk[l,k]
                end
            end
            j = i+k-1
            data[j] = profile[1:iprofile]
        end
    end
    return data
end


function loadmis3(ncvar::NCDatasets.Variable{T,2},
                  flag::NCDatasets.Variable{Tflag,2},
                  fillval,accepted_status_flag_values,
                  ncz::NCDatasets.Variable{Tz,2},
                  flag_z::NCDatasets.Variable{Tflagz,2},
                  fillval_z,accepted_status_flag_values_z) where {T,Tz,Tflag,Tflagz}
    nchunk = 10
    n_samples = size(ncvar,1)
    n_stations = size(ncvar,2)
#    n_stations = 100000
#    n_stations = 10000
    data = Vector{Vector{T}}(undef,n_stations);
    data_z = Vector{Vector{T}}(undef,n_stations);

    z = Vector{Vector{T}}(undef,n_stations);
    data_chunk = Array{T,2}(undef,(n_samples,nchunk))
    z_chunk = Array{Tz,2}(undef,(n_samples,nchunk))
    flag_chunk = Array{Tflag,2}(undef,(n_samples,nchunk))
    flag_z_chunk = Array{Tflag,2}(undef,(n_samples,nchunk))

    profile = Vector{T}(undef,n_samples)
    profile_z = Vector{T}(undef,n_samples)

    @inbounds for i = 1:nchunk:n_stations
        if (i-1) % (nchunk*1000) == 0
            println("$(i-1) out of $n_stations - $(100*(i-1)/n_stations) %")
        end
        nc = i:min(i+nchunk-1,n_stations)
        clen = length(nc)

        load!(ncvar,data_chunk,:,nc)
        load!(ncz,  z_chunk,   :,nc)
        load!(flag, flag_chunk,:,nc)
        load!(flag_z, flag_z_chunk,:,nc)

        for k = 1:clen
            iprofile = 0
            for l = 1:n_samples
                if ((data_chunk[l,k] != fillval)
                    && (z_chunk[l,k] != fillval_z)
                    && (flag_chunk[l,k] ∈ accepted_status_flag_values)
                    && (flag_z_chunk[l,k] ∈ accepted_status_flag_values_z)
                    )

                    iprofile = iprofile+1
                    profile[iprofile] = data_chunk[l,k]
                    profile_z[iprofile] = z_chunk[l,k]
                end
            end
            j = i+k-1
            data[j] = profile[1:iprofile]
            data_z[j] = profile_z[1:iprofile]
        end
    end
    return data,data_z
end





#data = @time loadmis(ncvar.var,fillval)
#data = @time loadmis2(ncvar.var,fillval)


#=

ll = 8529

@show extrema(length.(data))
@test data[ll][1:3] == [24.2607f0, 24.2608f0, 24.2525f0]


=#

function flatten_data(T,obsproflon,obsproflat,obsproftime,EDMO_CODE,LOCAL_CDI_ID,data,data_z)
    len = sum(length.(data))
    flat_lon = zeros(T,len)
    flat_lat = zeros(T,len)
    flat_time = Vector{DateTime}(undef,len)
    flat_ids = fill("",(len,))
    sel = trues(len)

    flat_data = Vector{T}(vcat(data...));
    flat_z = Vector{T}(vcat(data_z...));

    j = 0;

    for i = 1:length(data)
#    for i = 1:100
        jend = j+length(data[i]);
        obsid = "$(EDMO_CODE[i])-$(LOCAL_CDI_ID[i])"

        if (ismissing(obsproftime[i]) || ismissing(obsproflon[i])
            || ismissing(obsproflat[i]))
            sel[j+1:jend] .= false
        else
            flat_lon[j+1:jend] .= obsproflon[i]
            flat_lat[j+1:jend] .= obsproflon[i]
            flat_time[j+1:jend] .= obsproftime[i]
            flat_ids[j+1:jend] .= obsid
        end
        j = jend
    end

    return flat_data[sel],flat_lon[sel],flat_lat[sel],flat_z[sel],flat_time[sel],flat_ids[sel]
end

function flagvalues(attrib,accepted_status_flags)
    accepted_status_flags = ["good_value","probably_good_value"]

    flag_values = attrib["flag_values"]
    flag_meanings = attrib["flag_meanings"]::String
    if typeof(flag_meanings) <: AbstractString
        flag_meanings = split(flag_meanings)
    end

    accepted_status_flag_values = zeros(eltype(flag_values),length(accepted_status_flags))
    for i = 1:length(accepted_status_flags)
        tmp = findfirst(accepted_status_flags[i] .== flag_meanings)

        if tmp == nothing
            error("cannot recognise flag $(accepted_status_flags[i])")
        end
        accepted_status_flag_values[i] = flag_values[tmp]
    end

    return accepted_status_flag_values
end

"""
    obsvalue,obslon,obslat,obsdepth,obstime,obsids = load(T,fname,long_name;
         qv_flags = ["good_value","probably_good_value"])

Load all profiles in the file `fname` corresponding to netcdf variable with the
`long_name` attribute equal to long_name. `qv_flags` is a list of strings
with the quality flags to be kept. `obsids` is a vector of strings with the
EDMO code and local CDI id concatenated by a hypthen.
"""
function load(T,fname,long_name;
              qv_flags = ["good_value","probably_good_value"])

    accepted_status_flags = qv_flags
    ds = Dataset(fname);

    ncvar_LOCAL_CDI_ID = varbyattrib_first(ds,long_name = "LOCAL_CDI_ID")
    LOCAL_CDI_ID = DIVAnd.chararray2strings(ncvar_LOCAL_CDI_ID.var[:]);
    EDMO_CODE = varbyattrib_first(ds,long_name = "EDMO_CODE")[:];

    obsproflon = varbyattrib_first(ds,standard_name = "longitude")[:]
    obsproflat = varbyattrib_first(ds,standard_name = "latitude")[:]
    obsproftime = varbyattrib_first(ds,standard_name = "time")[:]

    ncvar = varbyattrib_first(ds,long_name = long_name);
    ncvar_z = varbyattrib_first(ds,long_name = "Depth");

    ncv_ancillary = NCDatasets.ancillaryvariables(ncvar,"status_flag").var
    ncv_ancillary_z = NCDatasets.ancillaryvariables(ncvar_z,"status_flag").var

    accepted_status_flag_values = flagvalues(ncv_ancillary.attrib,accepted_status_flags)
    accepted_status_flag_values_z = flagvalues(ncv_ancillary_z.attrib,accepted_status_flags)
    @show accepted_status_flag_values

    fillval = ncvar.attrib["_FillValue"]
    fillval_z = get(ncvar.attrib,"_FillValue",nothing)
    data,data_z = @time loadmis3(ncvar.var,ncv_ancillary,fillval,accepted_status_flag_values,
                                 ncvar_z.var,ncv_ancillary_z,fillval_z,accepted_status_flag_values_z)

    close(ds)

    obsvalue,obslon,obslat,obsdepth,obstime,obsids = flatten_data(T,obsproflon,obsproflat,obsproftime,EDMO_CODE,LOCAL_CDI_ID,data,data_z)
    return obsvalue,obslon,obslat,obsdepth,obstime,obsids
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

fname = "filename.nc"
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

obsvalue,obslon,obslat,obsdepth,obstime,obsids = load(T,fname,long_name; qv_flags = qv_flags)

@test obsvalue == [10.0, 20.0, 21.0, 12.0]
@test obslon == [1.0, 1.0, 2.0, 3.0]
@test obslat == [1.0, 1.0, 2.0, 3.0]
@test obstime == [DateTime(2000,1,1),DateTime(2000,1,1),DateTime(2000,1,2),DateTime(2000,1,3)]
@test obsdepth == [0.0, 1.0, 1.0, 0.0]
@test obsids == ["111-ax", "111-ax", "222-by", "333-cz"]
