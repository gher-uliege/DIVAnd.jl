
function deepest(array::AbstractArray{Union{Missing,T},3}) where T
    sz = size(array)
    array2d = Array{Union{Missing,T},3-1}(undef,sz[1],sz[2])
    array2d .= missing

    for k = 1:sz[3]
        for j = 1:sz[2]
            for i = 1:sz[1]
                if !ismissing(array[i,j,k])
                    # possibly overwrite previous value
                    array2d[i,j] = array[i,j,k]
                end
            end
        end
    end
    return array2d
end

function ncdeepest(nc,varname,suffix,thresholds_value)
    T = Float32
    ncvar = nc[varname * suffix]
    sz = size(ncvar)

    @assert length(sz) == 4
    dims = ("lon", "lat", "time")

    newvarname = varname * "_deepest" * suffix

    longname = ncvar.attrib["long_name"]
    valex = ncvar.attrib["_FillValue"]
    longname_deepest = "Deepest values of $(longname)"

    if thresholds_value !== nothing
        longname_deepest = longname_deepest * " masked using relative error threshold $(thresholds_value)"
    end

    # value of deepest layer
    @info("Creating new variable $(newvarname)")
    ncvar_deepest = defVar(nc, newvarname, T, dims,  attrib = OrderedDict(
        "_FillValue" => T(valex),
        "missing_value" => T(valex),
        "units" => ncvar.attrib["units"],
        "long_name" => "Deepest values of $(longname)",
    ))

    for n = 1:sz[4]
        field3D = nc[varname][:,:,:,n]
        ncvar_deepest[:,:,n] = deepest(field3D)
    end
end

function derived(datafile,varname;
                 error_thresholds = [("L1", 0.3), ("L2", 0.5)],
                 )

    NCDatasets.Dataset(datafile, "a") do nc
        ncvar = nc[varname]
        ncvar_relerr = nc[varname * "_relerr"]
        longname = ncvar.attrib["long_name"]
        T = Float32
        sz = size(ncvar)

        ncattrib = OrderedDict(ncvar.attrib)
        ncdeepest(nc,varname,"",nothing)
        dims = ("lon", "lat", "time")

        @show ncvar[10,3,:,1]
        @show nc[varname * "_deepest"][10,3,1]

        for (thresholds_name, thresholds_value) in error_thresholds
            newvarname = "$(varname)_deepest_$(thresholds_name)"
            ncattrib["long_name"] = "Deepest values of $(longname) masked using relative error threshold $(thresholds_value)"

            ncvar_Lx = defVar(nc, newvarname, T, dims,  attrib = ncattrib)
            for n = 1:sz[4]
                field = deepest(ncvar[:,:,:,n])
                relerr = deepest(ncvar_relerr[:,:,:,n])

                for ind in eachindex(relerr)
                    if !ismissing(relerr[ind])
                        if relerr[ind] .> thresholds_value
                            field[ind] = missing
                        end
                    end
                end

                ncvar_Lx[:,:,n] = field
            end
        end
        # add depth of deepest layer
        # Load the depth variable
        depth = nc["depth"][:]
        @info("Working on $(length(depth)) depth layers")

        # Load the 4D variable (we can load the first time instance,
        # as depth don't change with time)
        field3D = nc[varname][:,:,:,1]
        @info("size = $(size(field3D))")

        # Get missing value from field
        valex = nc[varname].attrib["missing_value"]

        newvarname = varname * "_deepest_depth"
        @info("Creating new variable $(newvarname)")
        ncvardeepestdepth = defVar(nc, newvarname, T, ("lon", "lat"))
        ncvardeepestdepth.attrib["long_name"] = "Deepest depth for $(varname)"
        ncvardeepestdepth.attrib["_FillValue"] = T(valex)
        ncvardeepestdepth.attrib["missing_value"] = T(valex)
        ncvardeepestdepth.attrib["units"] = "meters"
        ncvardeepestdepth.attrib["positive"] = "down"

      # Loop on depth: start from surface and go to the bottom
        # (I also add "abs" in case depth are negative, but not probable)
        depthindex = sortperm(abs.(depth))
        for idepth in depthindex
            @info("Depth index: $(idepth)")
            # Look for non-missing values at the considered depth
            nonmissing = (.!ismissing.(field3D[:,:,idepth]))
            @info("Found $(sum(nonmissing)) non missing values for depth $(depth[idepth])")

            # Write the variable
            ncvardeepestdepth[nonmissing] .= depth[idepth]
        end

        @info("Written new variable deepest depth")
    end
end



