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


function array2strings(obsids::Matrix{UInt8})
    obsid = Vector{String}(undef, size(obsids, 2))

    for i = 1:size(obsids, 2)
        id = view(obsids, :, i)
        index = findfirst(c -> c == 0, id)

        hasnonull = index == nothing

        obsid[i] = if hasnonull
            String(Char.(id))
        else
            String(Char.(view(id, 1:index-1)))
        end
    end

    return obsid
end

function loadobsid(ds, varname = "obsid")
    #return chararray2strings(ds[varname].var[:,:])
    data = Array{UInt8,2}(undef, size(ds["obsid"]))
    NCDatasets.load!(ds["obsid"].var, data, :, :)
    return array2strings(data)
end


"""
    obsvalue,obslon,obslat,obsdepth,obstime,obsid = loadobs(T,filename,varname)

Load the variable `varname` from the netCDF file `filename`.
Coordinates (the netCDF variables "obslon", "obslat", "obsdepth"),
time ("obstime") and identifiers ("obsids") will also be loaded.
Numeric output arguments will have the type `T`.

"""
function loadobs(T, filename, varname)

    Dataset(filename, "r") do ds
        time = nomissing(ds["obstime"][:])::Vector{DateTime}

        lon = Vector{T}(nomissing(ds["obslon"][:], NaN))
        lat = Vector{T}(nomissing(ds["obslat"][:], NaN))
        depth = Vector{T}(nomissing(ds["obsdepth"][:], NaN))
        value = Vector{T}(nomissing(ds[varname][:], NaN))

        obsid = loadobsid(ds, "obsid")

        return value, lon, lat, depth, time, obsid
    end
end
