"""
    value,lon,lat,depth,time,obsid = loadbigfile(filename)

Load data from the text file `filename` and returns vectors with
the value, longitude, latitude, depth and time (as DateTime).
A list string identifiers is also returned.
"""

function loadbigfile(fname)

    info("Loading data from 'big file' $(fname)")
    data = readlines(open(fname,"r"))
    nobs = length(data)

    lon = zeros(nobs)
    lat = zeros(nobs)
    depth = zeros(nobs)
    timeval = Array{DateTime}(nobs)
    value = zeros(nobs)
    id = Array{String}(nobs)

		mydate(x) = try
            DateTime(x)
           catch
		    @show x
            DateTime(Date(1900))
        end

    for i in 1:nobs
        rec = split(data[i])
        lon[i] = parse(Float64,rec[1])
        lat[i] = parse(Float64,rec[2])
        value[i] = parse(Float64,rec[3])
        depth[i] = parse(Float64,rec[4])



        #timeval[i] = DateTime(rec[10])
		timeval[i]=mydate(rec[10])
        id[i] = rec[11]
    end

    return value,lon,lat,depth,timeval,id
end


"""
    value,lon,lat,depth,time,obsid = loadobs(T,filename,varname)

Load the variable `varname` from the NetCDF file `filename`.
Coordinates (the NetCDF variables "obslon", "obslat", "obsdepth"),
time ("obstime") and identifies ("obsids") will also be loaded.
Numeric output arguments will have the type `T`.

"""
function loadobs(T,filename,varname)
    @inline function missingasNaN(v)
        v2 = fill(T(NaN),size(v))
        v2[.!ismissing.(v)] = v[.!ismissing.(v)]
        return v2
    end


    ds = Dataset(filename,"r")
    time = ds["obstime"][:].data;

    lon = missingasNaN(ds["obslon"][:])
    lat = missingasNaN(ds["obslat"][:])
    depth = missingasNaN(ds["obsdepth"][:])
    value = missingasNaN(ds[varname][:])


    obsids = ds["obsid"][:]

    obsid = Vector{String}(size(obsids,2))

    for i = 1:size(obsids,2)
        obsid[i] = strip(join(obsids[:,i]),'\0')
    end

    close(ds)
    return value,lon,lat,depth,time,obsid
end
