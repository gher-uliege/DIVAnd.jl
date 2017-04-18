function loadbigfile(fname)

    data = readlines(open(fname,"r"))
    nobs = length(data)

    lon = zeros(nobs)
    lat = zeros(nobs)
    depth = zeros(nobs)
    time = Array{DateTime}(nobs)
    value = zeros(nobs)
    id = Array{String}(nobs)


    for i in 1:nobs
        rec = split(data[i])
        lon[i] = parse(Float64,rec[1])
        lat[i] = parse(Float64,rec[2])
        value[i] = parse(Float64,rec[3])
        depth[i] = parse(Float64,rec[4])
        time[i] = DateTime(rec[10])
        id[i] = rec[11]
    end

    return value,lon,lat,depth,time,id
end
