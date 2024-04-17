# NetCDF file as produced by ODV

module NCODV

using Dates
using NCDatasets
using Missings
import ..chararray2strings


function varbyattrib_first(ds; kwargs...)
    vs = varbyattrib(ds; kwargs...)
    if length(vs) == 0
        str = join(["attribute '$k' equal to '$v'" for (k, v) in kwargs], " and ")
        error("No NetCDF variable found with $(str) in $(path(ds))")
    end

    if length(vs) > 1
        str = join(["attribute '$k' equal to '$v'" for (k, v) in kwargs], " and ")
        error("Several NetCDFs variable found with the $(str) in $(path(ds)). Loading this file is ambiguous. Please file an issue at https://github.com/gher-ulg/DIVAnd.jl/issues with the output of julia command `using NCDatasets; NCDataset(\"$(path(ds))\")` or the shell command `ncdump -h \"$(path(ds))\"` if this file has been produced by ODV.")
    end

    @debug begin
        str = join(["attribute '$k' equal to '$v'" for (k, v) in kwargs], " and ")
        @debug "use variable $(name(vs[1])) ($str)"
    end
    return vs[1]
end


_promote_Float64_or_more(x::Float32) = Float64(x)
_promote_Float64_or_more(x) = x

function decode_odv_years(data,fillvalue)
    t = similar(data, Union{DateTime,Missing})
    @inbounds for i in eachindex(data)
        if data[i] == fillvalue
            t[i] = missing
        else
            data_float64 = _promote_Float64_or_more(data[i])
            year = floor(Int,data_float64)
            yearlen = (Dates.isleapyear(year) ? 366 : 365)

            doy_ms = round(Int64,1000*24*60*60 * yearlen * (data_float64 - year))
            t[i] = DateTime(year,1,1) + Dates.Millisecond(doy_ms)
        end
    end
    return t
end


function alloc_chunk(flag,sz)
    if flag !== nothing
        Array{eltype(flag),2}(undef, sz)
    else
        nothing
    end
end

# # files always have variables with the long_name  "LOCAL_CDI_ID" and "EDMO_CODE" (all upper-case)
# # long_name for the primary variable to analysis are always P35 names
# # longitude, latitude and time (including dates) have the standard attribute "longitude", "latitude" and "time", respectively

function loadprof(
    ncvar::NCDatasets.Variable{T,2},
    flag,
    fillval,
    accepted_status_flag_values,

    ncz::NCDatasets.Variable{Tz,2},
    flag_z,
    fillval_z,
    accepted_status_flag_values_z,

    nctime,
    flag_time,
    fillval_time,
    accepted_status_flag_values_time;
    nchunk = 10
    ) where {T,Tz}

    @debug "accepted_status_flag_values_time: $accepted_status_flag_values_time, time name $(name(nctime))"

    n_samples = size(ncvar,1)
    n_stations = size(ncvar, 2)

    nchunk = min(nchunk,n_stations)
    #    n_stations = 100000
    #    n_stations = 10000
    data_chunk = Array{T,2}(undef, (n_samples, nchunk))
    flag_chunk = alloc_chunk(flag,(n_samples, nchunk))

    profile = Vector{T}(undef, n_samples)
    data = Vector{Vector{T}}(undef, n_stations)

    #z = Vector{Vector{T}}(undef, n_stations)
    z_chunk = Array{Tz,2}(undef, (n_samples, nchunk))
    flag_z_chunk = alloc_chunk(flag_z,(n_samples, nchunk))

    profile_z = Vector{T}(undef, n_samples)
    data_z = Vector{Vector{T}}(undef, n_stations)

    if nctime != nothing
        Ttime = eltype(nctime)
        time_chunk = Array{Ttime,2}(undef, (n_samples, nchunk))
        flag_time_chunk = alloc_chunk(flag_time,(n_samples, nchunk))

        profile_time = Vector{Float64}(undef, n_samples)
    end
    data_time = Vector{Vector{Float64}}(undef, n_stations)

    t0 = Base.time()
    @inbounds for i = 1:nchunk:n_stations
        t1 = Base.time()
        #if ((i - 1) % (nchunk * 1000) == 0) && (n_stations > 2 * nchunk)
        if t1 - t0 > 2
            println("$(i-1) out of $n_stations - $(100*(i-1)/n_stations) %")
            t0 = t1
        end
        nc = i:min(i + nchunk - 1, n_stations)
        clen = length(nc)

        NCDatasets.load!(ncvar, view(data_chunk,:,1:clen), :, nc)

        if flag != nothing
            NCDatasets.load!(flag, view(flag_chunk,:,1:clen), :, nc)
        end

        NCDatasets.load!(ncz, view(z_chunk,:,1:clen), :, nc)

        if flag_z != nothing
            NCDatasets.load!(flag_z, view(flag_z_chunk,:,1:clen), :, nc)
        end

        if nctime != nothing
            NCDatasets.load!(nctime, view(time_chunk,:,1:clen), :, nc)
            NCDatasets.load!(flag_time, view(flag_time_chunk,:,1:clen), :, nc)
        end

        for k = 1:clen
            iprofile = 0
            for l = 1:n_samples
                ok = true

                # check value
                ok = ok && (data_chunk[l, k] != fillval)

                if flag != nothing
                    ok = ok && (flag_chunk[l, k] ∈ accepted_status_flag_values)
                end

                # check depth
                ok = ok && (z_chunk[l, k] != fillval_z)

                if flag_z != nothing
                    ok = ok && (flag_z_chunk[l, k] ∈ accepted_status_flag_values_z)
                end

                # check time (for time series)
                if nctime != nothing
                    ok = ok && (time_chunk[l, k] != fillval_time)

                    if flag_time != nothing
                        ok = ok && (flag_time_chunk[l, k] ∈ accepted_status_flag_values_time)
                    end
                end

                if ok
                    iprofile = iprofile + 1
                    profile[iprofile] = data_chunk[l, k]
                    profile_z[iprofile] = z_chunk[l, k]

                    if nctime != nothing
                        profile_time[iprofile] = time_chunk[l, k]
                    end
                end
            end

            j = i + k - 1
            data[j] = profile[1:iprofile]
            data_z[j] = profile_z[1:iprofile]

            if nctime != nothing
                data_time[j] = profile_time[1:iprofile]
            end
        end
    end
    return data, data_z, data_time
end


function flatten_data(
    T,
    obsproflon,
    obsproflat,
    obsproftime,
    EDMO_CODE,
    LOCAL_CDI_ID,
    data,
    data_z,
)
    len = sum(length.(data))
    flat_lon = zeros(T, len)
    flat_lat = zeros(T, len)
    flat_ids = fill("", (len,))
    sel = trues(len)

    #flat_data = Vector{T}(vcat(data...))
    flat_data = Vector{T}(reduce(vcat,data))
    #flat_z = Vector{T}(vcat(data_z...))
    flat_z = Vector{T}(reduce(vcat,data_z))

    if ndims(obsproftime) == 2
        # time series
        flat_time = reduce(vcat,obsproftime)
    else
        # profile
        flat_time = Vector{DateTime}(undef, len)
    end

    j = 0

    for i = 1:length(data)
        #    for i = 1:100
        jend = j + length(data[i])
        obsid = "$(EDMO_CODE[i])-$(LOCAL_CDI_ID[i])"
        #@show obsproftime[i]

        notok = false
        notok = notok || ismissing(obsproflon[i])
        notok = notok || ismissing(obsproflat[i])

        if ndims(obsproftime) == 1
            ## check only for profiles
            notok = notok || ismissing(obsproftime[i])
        end

        if notok
            sel[j+1:jend] .= false
        else
            flat_lon[j+1:jend] .= obsproflon[i]
            flat_lat[j+1:jend] .= obsproflat[i]

            if ndims(obsproftime) == 1
                # profile
                flat_time[j+1:jend] .= obsproftime[i]
            end
            flat_ids[j+1:jend] .= obsid
        end
        j = jend
    end

    return flat_data[sel],
    flat_lon[sel],
    flat_lat[sel],
    flat_z[sel],
    flat_time[sel],
    flat_ids[sel]
end

function flagvalues(attrib, accepted_status_flags)
    flag_values = attrib["flag_values"]
    flag_meanings = attrib["flag_meanings"]::String
    if typeof(flag_meanings) <: AbstractString
        flag_meanings = split(flag_meanings)
    end

    accepted_status_flag_values = zeros(eltype(flag_values), length(accepted_status_flags))
    for i = 1:length(accepted_status_flags)
        tmp = findfirst(accepted_status_flags[i] .== flag_meanings)

        if tmp == nothing
            error("cannot recognise flag $(accepted_status_flags[i])")
        end
        accepted_status_flag_values[i] = flag_values[tmp]
    end

    return accepted_status_flag_values
end

function statusflags(ncvar,accepted_status_flags)
    tmp = NCDatasets.ancillaryvariables(ncvar, "status_flag")
    if tmp != nothing
        @debug "variable flag: $(name(tmp)) for $(name(ncvar))"
        accepted_status_flag_values =
            flagvalues(tmp.attrib, accepted_status_flags)
        @debug "accepted_status_flag_values: $accepted_status_flag_values"
        return tmp.var, accepted_status_flag_values
    else
        @warn "No ancillary variable status_flags found for variable $(name(ncvar)). All data points are assumed to be good."
        return nothing, nothing
    end
end

"""
    obsvalue,obslon,obslat,obsdepth,obstime,obsids = NCODV.load(T,fname,long_name;
         qv_flags = ["good_value","probably_good_value"],
         nchunk = 10
)

Load all profiles in the file `fname` corresponding to netCDF variable with the
`long_name` attribute equal to the parameter `long_name`. `qv_flags` is a list of strings
with the quality flags to be kept. The filtering of the quality flags is applied
to the data variables, time and depth coordinates. `obsids` is a vector of
strings with the EDMO code and local CDI id concatenated by a hyphen.

`nchunk` is the number of profiles read at a time. Large values of `nchunk` can increase
performance but requirer also more memory.

The variable with the following standard_name should exits:
* longitude
* latitude
* time

As well as the variable with the following long_name:
* LOCAL\\_CDI\\_ID
* EDMO\\_code or EDMO\\_CODE
* Depth

A guide how to export NetCDF files from ODV is available [here](https://github.com/gher-uliege/EMODnet-Chemistry-GriddedMaps/blob/main/doc/ODV_netCDF_export.md)

"""
function load(T, fname, long_name; qv_flags = ["good_value", "probably_good_value"],
              depth_name = "Depth",
              nchunk = 10)

    accepted_status_flags = qv_flags

    Dataset(fname) do ds
        nstations = Int(ds.dim["N_STATIONS"])
        LOCAL_CDI_ID = fill("",nstations)

        if length(varbyattrib(ds, long_name = "LOCAL_CDI_ID")) == 0
            @warn "No variable with the long_name attribute \'LOCAL_CDI_ID\' in $fname found. We use the empty string for LOCAL_CDI_ID instead."
        else
            ncvar_LOCAL_CDI_ID = varbyattrib_first(ds, long_name = "LOCAL_CDI_ID")

            if ndims(ncvar_LOCAL_CDI_ID) == 2
                LOCAL_CDI_ID = chararray2strings(Array(ncvar_LOCAL_CDI_ID.var))
            else
                @warn """The variable with the long_name attribute \'LOCAL_CDI_ID\' is expected to have two dimensions. For example the output of 'ncdump -h' of $fname should contain:
[...]
    char metavar4(N_STATIONS, STRING36) ;
                metavar4:long_name = "LOCAL_CDI_ID" ;
[...]

We use the empty string for LOCAL_CDI_ID instead.
    """
            end
        end

        EDMO_CODE = if length(varbyattrib(ds; long_name = "EDMO_code")) > 0
            Array(varbyattrib_first(ds, long_name = "EDMO_code"))
        else
            Array(varbyattrib_first(ds, long_name = "EDMO_CODE"))
        end

        obsproflon = Array(varbyattrib_first(ds, standard_name = "longitude"))
        obsproflat = Array(varbyattrib_first(ds, standard_name = "latitude"))

        # time for time series
        ncvar_time = nothing
        ncv_ancillary_time = nothing
        fillval_time = nothing
        accepted_status_flag_values_time = nothing
        vars_time_ISO8601 = varbyattrib(ds, long_name = "time_ISO8601")

        if length(vars_time_ISO8601) == 1
            # time series
            ncvar_time = vars_time_ISO8601[1]
            @assert ndims(ncvar_time) == 2
        else
            # profile
            obsproftime = Array(varbyattrib_first(ds, standard_name = "time"))
            @assert ndims(obsproftime) == 1
        end

        ncvar = varbyattrib_first(ds, long_name = long_name)
        ncvar_z = varbyattrib_first(ds, long_name = depth_name)

        @debug "variable: $(name(ncvar))"
        @debug "variable z: $(name(ncvar_z))"

        ncv_ancillary, accepted_status_flag_values = statusflags(ncvar,accepted_status_flags)
        ncv_ancillary_z, accepted_status_flag_values_z = statusflags(ncvar_z,accepted_status_flags)

        if ncvar_time !== nothing
            ncv_ancillary_time, accepted_status_flag_values_time = statusflags(ncvar_time,accepted_status_flags)
            fillval_time = get(ncvar_time.attrib, "_FillValue", nothing)
        end

        fillval = ncvar.attrib["_FillValue"]
        fillval_z = get(ncvar_z.attrib, "_FillValue", nothing)

        data, data_z, data_time = loadprof(
            ncvar.var,
            ncv_ancillary,
            fillval,
            accepted_status_flag_values,

            ncvar_z.var,
            ncv_ancillary_z,
            fillval_z,
            accepted_status_flag_values_z,

            #= these 4 variables are nothing for profiles =#
            (ncvar_time == nothing ? nothing : ncvar_time.var),
            ncv_ancillary_time,
            fillval_time,
            accepted_status_flag_values_time,
            nchunk = nchunk
        )


        if ncvar_time !== nothing
            time_units = ncvar_time.attrib["units"]
            @debug "time_units: $time_units (decode as fractional years)"
            @assert time_units == "years since 0000-01-01"
            obsproftime = decode_odv_years.(data_time,fillval_time)
        end

        obsvalue, obslon, obslat, obsdepth, obstime, obsids = flatten_data(
            T,
            obsproflon,
            obsproflat,
            obsproftime,
            EDMO_CODE,
            LOCAL_CDI_ID,
            data,
            data_z,
        )
        return obsvalue, obslon, obslat, obsdepth, obstime, obsids
    end
end

end
