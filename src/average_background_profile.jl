function background_profile_searchz(depthr; factor = 1)
    kmax = length(depthr)
    searchz = Vector{NTuple{2,eltype(depthr)}}(undef, kmax)

    for k = 1:kmax
        k0 = max(1, k - 1)
        k1 = min(kmax, k + 1)

        zcenter = depthr[k0] / 4 + depthr[k] / 2 + depthr[k1] / 4
        Δz = factor * (depthr[k1] - depthr[k0]) / 4
        searchz[k] = (zcenter - Δz, zcenter + Δz)
        # factor = 1, we have
        # searchz[k] = ((depthr[k0] + depthr[k])/2,(depthr[k1] + depthr[k])/2)
    end
    return searchz
end

function simple_background_profile(
    obsdepth,
    obsvalue,
    depthr;
    epsilon2 = ones(size(obsvalue)),
    searchz = background_profile_searchz(depthr),
)

    kmax = length(depthr)
    profile = zeros(kmax)
    w = 1 ./ epsilon2

    for k = 1:kmax
        z_min, z_max = searchz[k]

        sel = (z_min .<= obsdepth .<= z_max) .& (isfinite.(obsvalue))
        profile[k] = mean(w[sel] .* obsvalue[sel]) / mean(w[sel])
    end

    return profile
end


"""
    average_background_profile(
    background_filename, (lonr,latr,depthr,TS), (obslon, obslat, obsdepth, obstime), obsvalue,
    epsilon2,
    varname;
    transform = DIVAnd.Anam.notransform(),
    searchz = background_profile_searchz(depthr),

Compute the average background profile by averaging the observations
within a distance given by `searchz` and for each time instance defined in
the time selector `TS`.
"""
function average_background_profile(
    background_filename,
    grid_range,
    (obslon, obslat, obsdepth, obstime),
    obsvalue,
    epsilon2,
    varname;
    transform = DIVAnd.Anam.notransform(),
    searchz = background_profile_searchz(grid_range[3]), # depthr
    chunksizes = [100, 100, 1, 1],
    checksum = :fletcher32,
    deflatelevel = 5,
)

    (lonr, latr, depthr, TS) = grid_range
    # anamorphosis transform
    trans, invtrans = transform

    sz = (length(lonr), length(latr), length(depthr))

    ds = Dataset(background_filename, "c")

    ds.dim["lon"] = sz[1]
    ds.dim["lat"] = sz[2]
    ds.dim["depth"] = sz[3]
    ds.dim["time"] = length(TS)


    chunksizes[1:3] = min.(chunksizes[1:3], collect(sz))

    ncvar = defVar(
        ds,
        varname,
        Float32,
        ("lon", "lat", "depth", "time"),
        deflatelevel = deflatelevel,
        checksum = checksum,
        chunksizes = chunksizes,
    )


    for timeindex = 1:length(TS)
        @info "Time step $(timeindex) / $(length(TS))"
        # select observation to be used for the time instance timeindex
        sel = DIVAnd.select(TS, timeindex, obstime)


        # apply the transformation
        value_trans = trans.(obsvalue[sel])
        # some transformation (e.g. log) can produce -Inf, set these to NaN
        value_trans[.!isfinite.(value_trans)] .= NaN

        bp = simple_background_profile(
            obsdepth[sel],
            value_trans,
            depthr;
            epsilon2 = epsilon2[sel],
            searchz = searchz,
        )

        background = zeros(sz)
        background[:, :, :] .= reshape(bp, (1, 1, :))

        ncvar[:, :, :, timeindex] = background
    end

    close(ds)

    return nothing
end
