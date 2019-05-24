# NetCDF file as produced by ODV

module NCODV

using NCDatasets
using Missings
import ..chararray2strings

if VERSION >= v"0.7.0-beta.0"
    using Dates
else
    using Compat: @info, @warn, @debug
end
using Compat


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

function loadprof(ncvar::NCDatasets.Variable{T,2},
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
        if ((i-1) % (nchunk*1000) == 0) && (n_stations > 2*nchunk)
            println("$(i-1) out of $n_stations - $(100*(i-1)/n_stations) %")
        end
        nc = i:min(i+nchunk-1,n_stations)
        clen = length(nc)

        NCDatasets.load!(ncvar,data_chunk,:,nc)
        NCDatasets.load!(ncz,  z_chunk,   :,nc)
        NCDatasets.load!(flag, flag_chunk,:,nc)
        NCDatasets.load!(flag_z, flag_z_chunk,:,nc)

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
            flat_lat[j+1:jend] .= obsproflat[i]
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
    obsvalue,obslon,obslat,obsdepth,obstime,obsids = NCODV.load(T,fname,long_name;
         qv_flags = ["good_value","probably_good_value"])

Load all profiles in the file `fname` corresponding to netcdf variable with the
`long_name` attribute equal to the parameter `long_name`. `qv_flags` is a list of strings
with the quality flags to be kept. `obsids` is a vector of strings with the
EDMO code and local CDI id concatenated by a hyphen.
"""
function load(T,fname,long_name;
              qv_flags = ["good_value","probably_good_value"])

    accepted_status_flags = qv_flags

    Dataset(fname) do ds
        ncvar_LOCAL_CDI_ID = varbyattrib_first(ds,long_name = "LOCAL_CDI_ID")
        LOCAL_CDI_ID = chararray2strings(ncvar_LOCAL_CDI_ID.var[:]);

        EDMO_CODE =
            if length(varbyattrib(ds; long_name = "EDMO_code")) > 0
                varbyattrib_first(ds,long_name = "EDMO_code")[:]
            else
                varbyattrib_first(ds,long_name = "EDMO_CODE")[:]
            end

        obsproflon = varbyattrib_first(ds,standard_name = "longitude")[:]
        obsproflat = varbyattrib_first(ds,standard_name = "latitude")[:]
        obsproftime = varbyattrib_first(ds,standard_name = "time")[:]

        ncvar = varbyattrib_first(ds,long_name = long_name);
        ncvar_z = varbyattrib_first(ds,long_name = "Depth");

        ncv_ancillary = NCDatasets.ancillaryvariables(ncvar,"status_flag").var
        ncv_ancillary_z = NCDatasets.ancillaryvariables(ncvar_z,"status_flag").var

        accepted_status_flag_values = flagvalues(ncv_ancillary.attrib,accepted_status_flags)
        accepted_status_flag_values_z = flagvalues(ncv_ancillary_z.attrib,accepted_status_flags)
        @debug accepted_status_flag_values

        fillval = ncvar.attrib["_FillValue"]
        fillval_z = get(ncvar.attrib,"_FillValue",nothing)
        data,data_z = loadprof(ncvar.var,ncv_ancillary,fillval,accepted_status_flag_values,
                               ncvar_z.var,ncv_ancillary_z,fillval_z,accepted_status_flag_values_z)

        obsvalue,obslon,obslat,obsdepth,obstime,obsids = flatten_data(T,obsproflon,obsproflat,obsproftime,EDMO_CODE,LOCAL_CDI_ID,data,data_z)
        return obsvalue,obslon,obslat,obsdepth,obstime,obsids
    end
end

end
