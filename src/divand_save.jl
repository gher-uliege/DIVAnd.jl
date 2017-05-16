function divand_save(filename,mask,varname,fi)

    sz = size(mask)
    dims = [NcDim("longitude",sz[1]),
            NcDim("latitude",sz[2]),
            NcDim("depth",sz[3]),
            NcDim("time",sz[4])]

    @show filename
    nc = NetCDF.create(filename,NcVar(varname,dims))
    nc[varname][:,:,:,:] = fi
    NetCDF.close(nc)

    return nothing
end