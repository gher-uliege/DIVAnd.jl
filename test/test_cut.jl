using PyPlot
using Missings



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



function cut(filename2,filename_cut,varname,maskkeep)
    ds = NCDataset(filename2)
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

    ds_cut = NCDataset(filename_cut,"c",attrib = ds.attrib)

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

    #deflatelevel = 5
    #chunksizes = [100, 100, 1, 1]

    for (varname,ncvar) in ds
        #storage,chunksizes = chunking(ncvar)
        #checksum(ncvar)
        #isshuffled,isdeflated,deflatelevel = deflate(ncvar)
        #@show chunksizes,storage
        dn = dimnames(ncvar)

        kwargs = Dict()
        if length(dn) > 3
            #        kwargs[:chunksizes] = collect(min.(chunk.(dn),size(ncvar)))
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
    close(ds_cut)
    close(ds)
end

lonr = -0.2:0.1:1.1
latr = -0.2:0.1:1.1

polygon_lon = [1.05, 1.05, 0]
polygon_lat = [0, 1.05, 1.05]

sz = (length(lonr),length(latr))
field = lonr .+ latr'

ds = NCDataset(filename2)
lonr = ds["lon"][:]
latr = ds["lat"][:]
close(ds)

polygon_lon = [5., 10, 11]
polygon_lat = [41, 41, 43.]
maskkeep = inpolygon(polygon_lon,polygon_lat,lonr,latr)
filename_cut = "/tmp/cut.nc"

cut(filename2,filename_cut,varname,maskkeep)

#=
clf()
subplot(2,1,1);
#pcolormesh(lonr,latr,Float64.(maskkeep)'); colorbar()
pcolormesh(lonr,latr,field'); 
plot(polygon_lon[[1:end; 1]],polygon_lat[[1:end; 1]],"k-")
cl = extrema(field)
ax = axis()
clim(cl)
colorbar()

subplot(2,1,2);
pcolormesh(lonr[i],latr[j],nomissing(fieldkeep,NaN)');
plot(polygon_lon[[1:end; 1]],polygon_lat[[1:end; 1]],"k-")
clim(cl)
axis(ax)
colorbar()
=#
