function _inpolygon(X,(testx,testy))
    c = false
    nvert = size(X,1)
    j = nvert

    for i = 1:nvert
        if ( ((X[i,2]>testy) != (X[j,2]>testy)) &&
             (testx < (X[j,1]-X[i,1]) * (testy-X[i,2]) / (X[j,2]-X[i,2]) + X[i,1]) )
            c = !c
        end
        j = i
    end
    return c
end

function inpolygon(polygon_lon,polygon_lat,lonr,latr)
    sz = (length(lonr),length(latr))
    maskkeep = trues(sz)

    X = hcat(polygon_lon,polygon_lat)

    for j = 1:sz[2]
        for i = 1:sz[1]
            maskkeep[i,j] = _inpolygon(X,(lonr[i],latr[j]))
        end
    end
    return maskkeep
end


function _crop(maskkeep,field)
    fieldsub = allowmissing(field[i,j])
    fieldsub[.!maskkeepsub] .= missing
    return fieldsub
end



function maskout!(dn,maskkeepsub,field)
    if length(dn) >= 2
        if dn[1:2] == ("lon","lat")
            for ind in CartesianIndices((1:size(field,3), 1:size(field,4)))
                tmp = @view field[:,:,ind]
                tmp[.!maskkeepsub] .= missing
            end
        end
    end
end

"""
    DIVAnd.cut(filename,varname,filename_cut,polygon_lon,polygon_lat)

Exclude all grid cells points of the netcdf file `filename` with the
variable `varname` to the grid cell included in the polygon
whose vertices are defined by the vector of longitude values `polygon_lon`
and the vector of latitude values `polygon_lat`. Values outside of this polygone
will be croped (if possible) or masked (if croping is not possible).

Only the NetCDF variable with the dimensions `lon` and `lat` will be
modified. All other variarables (in particular `obslon`, `obslat`, ...)
will not be changed.

Example

```julia
filename2 = "Water_body_dissolved_oxygen_concentration_monthly.nc"
polygon_lon = [-20., 20, 20, -19]
polygon_lat = [21, 21, 53., 52.]
varname = "Water body dissolved oxygen concentration"
filename_cut = "cut.nc"
cut(filename2,varname,filename_cut,polygon_lon,polygon_lat)
```

"""
function cut(filename2,varname,filename_cut,
             polygon_lon::AbstractVector,
             polygon_lat::AbstractVector; kwargs...)

    (lonr,latr) =
        NCDataset(filename2) do ds
            ds["lon"][:], ds["lat"][:]
        end

    maskkeep = inpolygon(polygon_lon,polygon_lat,lonr,latr)
    @info "keeping: $(count(maskkeep))"
    cut(filename2,varname,filename_cut,maskkeep; kwargs...)
end

function cut(filename2,varname,filename_cut,maskkeep::AbstractArray{Bool,2};
             compress = true
             )
    NCDataset(filename2) do ds
        lonr = ds["lon"][:]
        latr = ds["lat"][:]

        maskkeep1 = any(maskkeep,dims=2)[:,1]
        i = findfirst(maskkeep1):findlast(maskkeep1)

        maskkeep2 = any(maskkeep,dims=1)[1,:]
        j = findfirst(maskkeep2):findlast(maskkeep2)

        maskkeepsub = maskkeep[i,j]

        #fieldkeep = _crop(maskkeep,field)

        if isfile(filename_cut)
            rm(filename_cut)
        end

        NCDataset(filename_cut,"c",attrib = ds.attrib) do ds_cut

            # set dimenions
            for (dimname,dimlen) in ds.dim
                ds_cut.dim[dimname] =
                    if dimname == "lon"
                        length(i)
                    elseif dimname == "lat"
                        length(j)
                    else
                        dimlen
                    end
            end

            function slice_(dimname)
                if dimname == "lon"
                    i
                elseif dimname == "lat"
                    j
                else
                    Colon()
                end
            end


            function chunk(dimname)
                if dimname == "lon"
                    100
                elseif dimname == "lat"
                    100
                else
                    1
                end
            end

            for (varname,ncvar) in ds
                dn = dimnames(ncvar)
                @debug "varname: $varname"
                @debug "compress: $compress"
                kwargs = Dict()

                if compress
                    storage,chunksizes = chunking(ncvar)
                    isshuffled,isdeflated,deflatelevel = deflate(ncvar)
                    sizecut = ntuple(i -> ds_cut.dim[dn[i]],length(dn))

                    @debug "sizecut: $(sizecut)"
                    chunksizes  = min.(chunksizes,sizecut)

                    @debug "checksum: $(checksum(ncvar))"
                    @debug "chunksizes: $chunksizes"
                    @debug "storage: $chunksizes"
                    @debug "deflatelevel: $deflatelevel"
                    @debug "size(ncvar): $(size(ncvar))"


                    if !(varname in ["lon","lat","time","depth","climatology_bounds"])
                        kwargs[:checksum] = checksum(ncvar)
                        if storage == :chunked
                            kwargs[:chunksizes] = chunksizes
                        end
                        if isdeflated
                            kwargs[:deflatelevel] = deflatelevel
                        end
                    end
                    @debug "arguments: $(kwargs)"
                end

                #@show chunksizes,size(ncvar)

                T = eltype(ncvar.var)
                ncvar_cut = defVar(
                    ds_cut,varname,T,dn;
                    kwargs...,
                    #        deflatelevel = deflatelevel,
                    #        chunksizes = chunksizes,
                    #        checksum = checksum(ncvar),
                    attrib = OrderedDict(ncvar.attrib))

                indices = slice_.(dn)
                indices_write = ntuple(i -> Colon(), length(dn))

                @info "Slicing $varname[$indices]"

                if "time" in dn
                    for n = 1:ds.dim["time"]
                        data = ncvar[indices[1:end-1]...,n]
                        maskout!(dn,maskkeepsub,data)
                        ncvar_cut[indices_write[1:end-1]...,n] = data
                    end
                else
                    data = ncvar[indices...]
                    maskout!(dn,maskkeepsub,data)
                    ncvar_cut[indices_write...] = data
                end
            end
        end
    end
end
