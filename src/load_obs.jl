"""
    value,lon,lat,depth,time,obsid = loadbigfile(filename)

Load data from the text file `filename` and returns vectors with
the value, longitude, latitude, depth and time (as DateTime).
A list string identifiers is also returned.
"""
function loadbigfile(fname)

    @info "Loading data from 'big file' $(fname)"
    data = readlines(open(fname, "r"))
    nobs = length(data)

    lon = zeros(nobs)
    lat = zeros(nobs)
    depth = zeros(nobs)
    timeval = Array{DateTime}(undef, nobs)
    value = zeros(nobs)
    id = Array{String}(undef, nobs)

    mydate(x) =
        try
            DateTime(x)
        catch
            @show x
            DateTime(Date(1900))
        end

    for i = 1:nobs
        rec = split(data[i])
        lon[i] = parse(Float64, rec[1])
        lat[i] = parse(Float64, rec[2])
        value[i] = parse(Float64, rec[3])
        depth[i] = parse(Float64, rec[4])



        #timeval[i] = DateTime(rec[10])
        timeval[i] = mydate(rec[10])
        id[i] = rec[11]
    end

    return value, lon, lat, depth, timeval, id
end

function loadobsid(filename::AbstractString, varname = "obsid")
    Dataset(filename, "r") do ds
        return loadobsid(ds, varname)
    end
end


function chararray2strings(obsids::Matrix{Char})
    obsid = Vector{String}(undef, size(obsids, 2))

    for i = 1:size(obsids, 2)
        id = view(obsids, :, i)
        index = findfirst(c -> c == '\0', id)

        hasnonull = index == nothing

        obsid[i] = if hasnonull
            String(id)
        else
            String(view(id, 1:index-1))
        end
    end

    return obsid
end


function array2strings!(obsids::AbstractMatrix{UInt8},obsid)
    len,nobs = size(obsids)

    previous_id = zeros(UInt8,len)
    previous_obsid = ""

    for i = 1:nobs
        id = view(obsids, :, i)

        if (id == previous_id) && (i > 1)
            # reuse previous ids to save memory
            obsid[i] = previous_obsid
        else
            # create a new string
            index = findfirst(c -> c == 0, id)

            hasnonull = index == nothing

            previous_obsid = if hasnonull
                String(Char.(id))
            else
                String(Char.(view(id, 1:index-1)))
            end
            obsid[i] = previous_obsid
            previous_id .= id
        end

    end

    return obsid
end

function array2strings(obsids::AbstractMatrix{UInt8})
    obsid = Vector{String}(undef, size(obsids, 2))
    return array2strings!(obsids,obsid)
end


function loadobsid(ds, varname = "obsid"; chunksize = 1_000_000)
    len,nobs = size(ds["obsid"]) :: Tuple{Int,Int}
    obsids = Vector{String}(undef,nobs)

    data = Array{UInt8,2}(undef, (len,chunksize))

    for j = 1:chunksize:nobs
        i = j:min(j+chunksize-1,nobs)
        k = i .- (j-1)

        data .= 0
        NCDatasets.load!(ds["obsid"].var, view(data,:,k), :, i)
        array2strings!(view(data,:,k),view(obsids,i))
    end
    return obsids
end


"""
    obsvalue,obslon,obslat,obsdepth,obstime,obsid = loadobs(T,filename,varname)

Load the variable `varname` from the netCDF file `filename`.
Coordinates (the netCDF variables "obslon", "obslat", "obsdepth"),
time ("obstime") and identifiers ("obsids") will also be loaded.
Numeric output arguments will have the type `T`.

"""
function loadobs(T, filename, varname; chunksize = 1_000_000)

    Dataset(filename, "r") do ds
        time = nomissing(ds["obstime"][:])::Vector{DateTime}

        lon = Vector{T}(nomissing(ds["obslon"][:], NaN))
        lat = Vector{T}(nomissing(ds["obslat"][:], NaN))
        depth = Vector{T}(nomissing(ds["obsdepth"][:], NaN))
        value = Vector{T}(nomissing(ds[varname][:], NaN))

        obsid = loadobsid(ds, "obsid"; chunksize = chunksize)

        return value, lon, lat, depth, time, obsid
    end
end
