
"""
    bx,by,b = DIVAnd.extract_bath(bath_name,isglobal,xi,yi)

Extract the bathymetry from the NetCDF file `bathname`. The parameter `isglobal`
is true if the NetCDF file covers the whole globe and thus the last longitude
point can be considered to be right next to the first longitude point.
`xi` and `yi` are vectors defining the bounding box of the data.
No interpolation is performed.

**Convention:** b is positive in the water and negative in the air.

The NetCDF file is expected to have the one dimensional variables `lon` and `lat` with the
longitude (degrees East) and latitude (degrees North) and the two
dimentional array `bat` with the digital terrain model
(negative in water and positive above water). The order of the dimension
should follow be: longitude and then latitude in
[Column-major ordering](https://en.wikipedia.org/wiki/Row-_and_column-major_order)
(or latitude and then longitude if the tool `ncdump` is used,
which is based on Row-major ordering).

Example of the output of `ncdump -h`:
```
netcdf gebco_30sec_8 {
dimensions:
	lat = 2702 ;
	lon = 5400 ;
variables:
	double lat(lat) ;
		lat:long_name = "Latitude" ;
		lat:standard_name = "latitude" ;
		lat:units = "degrees_north" ;
	double lon(lon) ;
		lon:long_name = "Longitude" ;
		lon:standard_name = "longitude" ;
		lon:units = "degrees_east" ;
	float bat(lat, lon) ;
		bat:long_name = "elevation above sea level" ;
		bat:standard_name = "height" ;
		bat:units = "meters" ;

// global attributes:
		:title = "GEBCO" ;
}
```
"""
function extract_bath(bath_name, isglobal, xi, yi)

    #@info "Extracting bathymetry from file: $(bath_name)"
    # if isglobal == true
    #@info "Working with a global bathymetry";
    # end;

    # TODO: make it work for 2-element vector with the bounding box
    dxi = xi[2] - xi[1]
    dyi = yi[2] - yi[1]

    Dataset(bath_name) do nc
        x = nc["lon"].var[:]
        y = nc["lat"].var[:]

        if (maximum(yi) > maximum(y)) || (minimum(yi) < minimum(y))
            error("data sets '$bath_name' covers only latitude range $(extrema(y)), while latitude range $(extrema(yi)) is requested")
        end
        rx = x[2] - x[1]
        ry = y[2] - y[1]
        X0 = x[1]
        Y0 = y[1]

        redx = dxi / rx
        redy = dyi / ry

        # do not range check i0 and i1 now because of wrapping when bathymetry
        # is global
        i0 = floor(Int,(xi[1]-dxi-X0)/rx) + 1
        i1 = ceil(Int,(xi[end]+dxi-X0)/rx) + 1
        i = i0:i1

        j0 = max(floor(Int,(yi[1]-dyi-Y0)/ry)+1,1);
        j1 = min(ceil(Int,(yi[end]+dyi-Y0)/ry)+1,length(y));
        j = j0:j1
        b = zeros(length(i), length(j))

        if isglobal
            i2 = mod.(i .- 1, size(nc["bat"], 1)) .+ 1
            jumps = [0; findall(abs.(i2[2:end] - i2[1:end-1]) .> 1); length(i2)]

            for l = 1:length(jumps)-1
                b[(jumps[l]+1):jumps[l+1], :] = nc["bat"].var[
                    i2[jumps[l]+1]:i2[jumps[l+1]],
                    j,
                ]
            end
        else
            i = max(minimum(i),1):min(maximum(i),length(x))
            j = max(minimum(j),1):min(maximum(j),length(y))
            b[:, :] = nc["bat"].var[i, j]
        end

        bx = X0 .+ rx * (i .- 1)
        by = Y0 .+ ry * (j .- 1)

        return bx, by, -b
    end
end

"""
    xi,yi,bath = DIVAnd.load_bath(bath_name,isglobal,xi,yi)

Load the bathymetry from the netCDF file `bathname`. The parameter `isglobal` is true if the NetCDF file covers
the whole globe and thus the last longitude point can be considered to be right next to the first longitude point.
`xi` and `yi` are vectors containing the longitude and latitude grid onto which the bathymetry should be
interpolated.

"""
function load_bath(bath_name, isglobal, xi, yi)

    bx, by, b = extract_bath(bath_name, isglobal, xi, yi)

    # hack
    #b = -b
    #b[isnan.(b)] = 10
    # end hack
    #@show extrema(b)


    # # convolution
    # Fx = round(redx);
    # Fy = round(redy);

    # Fi,Fj = ndgrid(-Fx:Fx,-Fy:Fy);
    # F = exp(-2* ((Fi/redx).^2 + (Fj/redy).^2));
    # F = F/sum(F[:]);

    # # mask as float

    # maskf = Int.(mask) + 0.
    # #m = conv2(double(mask),F,"same");
    # m = conv2(maskf,F);
    Xi, Yi = ndgrid(xi, yi)

    itp = interpolate((bx, by), b, Gridded(Linear()))
    bi = itp(xi, yi)

    return xi, yi, bi
end



"""
    xi,yi,mask = load_mask(bath_name,isglobal,xi,yi,level::Number)

Generate a land-sea mask based on the topography from the NetCDF file
`bathname`. The parameter `isglobal` is true if the NetCDF file covers the whole globe and
thus the last longitude point can be considered to be right next to the first longitude point.
`xi` and `yi` are vectors containing the longitude and latitude grid onto which the bathymetry should be
interpolated.

**Convention:** in the water, `level` is positive and in the air `level` is negative.

"""
function load_mask(bath_name, isglobal, xi, yi, level::Number)

    #@info "Creating land-sea mask on level: $(level)"

    bx, by, b = extract_bath(bath_name, isglobal, xi, yi)

    # hack
    #b = -b
    #b[isnan.(b)] = 10
    # end hack
    #@show extrema(b)
    mask = b .> level


    # # convolution
    # Fx = round(redx);
    # Fy = round(redy);

    # Fi,Fj = ndgrid(-Fx:Fx,-Fy:Fy);
    # F = exp(-2* ((Fi/redx).^2 + (Fj/redy).^2));
    # F = F/sum(F[:]);

    # # mask as float

    # maskf = Int.(mask) + 0.
    # #m = conv2(double(mask),F,"same");
    # m = conv2(maskf,F);
    Xi, Yi = ndgrid(xi, yi)

    itp = interpolate((bx, by), Int.(mask), Gridded(Linear()))
    mif = itp(xi, yi)

    mi = mif .> 1 / 2

    return xi, yi, mi
end


function load_mask(bath_name, isglobal, x, y, levels::AbstractVector)
    data = [load_mask(bath_name, isglobal, x, y, level) for level in levels]

    catdata = cat([d[3] for d in data]..., dims = 3)
    return data[1][1], data[1][2], catdata
end
