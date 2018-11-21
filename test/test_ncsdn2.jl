using NCDatasets
using Missings
import DIVAnd
if VERSION >= v"0.7.0-beta.0"
    using Test
else
    using Base.Test
end

NCSDN = DIVAnd.NCSDN
#=
basedir = joinpath(@__DIR__,"..","data")

fname = joinpath(basedir,"netCDF_vertical_profiles_ctd.nc")

param = "TEMPPR01"

T = Float64
fnames = [fname]

value,lon,lat,z,time,ids = NCSDN.load(T,fnames,param;
     qualityflags = [NCSDN.GOOD_VALUE, NCSDN.PROBABLY_GOOD_VALUE])

fnames = [
    joinpath(basedir,"netCDF_timeseries_tidegauge.nc"),
    joinpath(basedir,"netCDF_timeseries_tidegauge_with_instrument.nc"),
    joinpath(basedir,"netCDF_trajectory_meteorological_data.nc"),
    joinpath(basedir,"netCDF_trajectory_tsg_with_instrument.nc"),
    joinpath(basedir,"netCDF_vertical_profiles_ctd.nc"),
    joinpath(basedir,"netCDF_vertical_profiles_ctd_with_instruments.nc"),
    joinpath(basedir,"netCDF_vertical_profiles_xbt_with_fall_rate_and_instruments.nc")]

value,lon,lat,z,time,ids = NCSDN.load(T,fnames,param;
     qualityflags = [NCSDN.GOOD_VALUE, NCSDN.PROBABLY_GOOD_VALUE])


param = "TEMPET01"

value,lon,lat,z,time,ids = NCSDN.load(T,fnames,param;
    qualityflags = [NCSDN.GOOD_VALUE, NCSDN.PROBABLY_GOOD_VALUE])

=#




fnames = ["/home/ulg/gher/abarth/Downloads/data_from_SDN_2017-11_TS_profiles_non-restricted_med.nc"]
fname = fnames[1]
accepted_status_flags = "good_value","probably_good_value";


long_name = "ITS-90 water temperature"

ds = Dataset(fname);


function varbyattrib_first(ds; kwargs...)
    vs = varbyattrib(ds; kwargs...)
    if length(vs) == 0
        str = join(["attribute '$k' equal to '$v'" for (k,v) in kwargs]," and ")
        error("No NetCDF variable found with $(str) in $(path(ds))")
    end
    return vs[1]
end

# LOCAL_CDI_ID = varbyattrib_first(ds,long_name = "LOCAL_CDI_ID");
# EDMO_CODE = varbyattrib_first(ds,long_name = "EDMO_CODE");

# obslat = varbyattrib_first(ds,standard_name = "latitude")[:]
# obslon = varbyattrib_first(ds,standard_name = "longitude")[:]
# obstime = varbyattrib_first(ds,standard_name = "time")[:]


# @show extrema(skipmissing(obstime))

# # files always hava variable with the long_name  "LOCAL_CDI_ID" and "EDMO_CODE" (all upper-case)
# # long_name for the primary variable to analysis are always P35 names
# # longitude, latitude and time (including dates) have the standard attribute "longitude", "latitude" and "time" respectively

ncvar = varbyattrib_first(ds,long_name = long_name);
ncvar_z = varbyattrib_first(ds,long_name = "Depth");

# n_stations = ds.dim["N_STATIONS"]
# T = Float32
# data = Vector{Vector{T}}(undef,n_stations);

# for i = 1:10
#     prof = NCDatasets.filter(ncvar, :,i, accepted_status_flags = accepted_status_flags);

#     data[i] = collect(skipmissing(prof))
#     #append!(prof,skipmissing(prof))
# end


@inline function load!(ncvar::NCDatasets.Variable{T,2}, data, i::Colon,j::UnitRange) where T
    n_samples = size(ncvar,1)
    # reversed and 0-based
    NCDatasets.nc_get_vara(ncvar.ncid,ncvar.varid,[first(j)-1,0],[length(j),n_samples],data)
#    @show datac[1:10],size(datac)
end


@inline function load!(ncvar::NCDatasets.Variable{T,2}, data, indices::Union{Integer, UnitRange, Colon}...) where T
    ind = to_indices(ncvar,indices)
    start = [first(ind[2])-1,first(ind[1])-1]
    count = [length(ind[2]),length(ind[1])]
    stride = [step(ind[2]),step(ind[1])]
    NCDatasets.nc_get_vars(ncvar.ncid,ncvar.varid,start,count,stride,data)
end

function loadmis(ncvar::NCDatasets.Variable{T,2},fillval) where T
    n_samples,n_stations = size(ncvar)
    n_stations = 100000
    data = Vector{Vector{T}}(undef,n_stations);
    tmp = Array{T,2}(undef,(n_samples,1))
    
    @inbounds for i = 1:n_stations
        # reversed and 0-based
        NCDatasets.nc_get_vars(ncvar.ncid,ncvar.varid,[i-1,0],[1,n_samples],[1,1],tmp)
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



function loadmis3(ncvar::NCDatasets.Variable{T,2},fillval,ncz::NCDatasets.Variable{Tz,2},fillval_z,flag::NCDatasets.Variable{Tflag,2},accepted_status_flag_values) where {T,Tz,Tflag}
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

        for k = 1:clen
            iprofile = 0
            for l = 1:n_samples
                if ((data_chunk[l,k] != fillval)
                    && (z_chunk[l,k] != fillval_z)
                    && (flag_chunk[l,k] ∈ accepted_status_flag_values)
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




ncv_ancillary = NCDatasets.ancillaryvariables(ncvar,"status_flag").var

accepted_status_flags = ["good_value","probably_good_value"]

flag_values = ncv_ancillary.attrib["flag_values"]
flag_meanings = ncv_ancillary.attrib["flag_meanings"]::String
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

@show accepted_status_flag_values

fillval = ncvar.attrib["_FillValue"]
fillval_z = get(ncvar.attrib,"_FillValue",nothing)


#data = @time loadmis(ncvar.var,fillval)
#data = @time loadmis2(ncvar.var,fillval)


data,data_z = @time loadmis3(ncvar.var,fillval,ncvar_z.var,fillval_z,ncv_ancillary,accepted_status_flag_values)

ll = 8529 

@show extrema(length.(data))
@test data[ll][1:3] == [24.2607f0, 24.2608f0, 24.2525f0]