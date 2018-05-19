
"""
    bx,by,b = divand.extract_bath(bath_name,isglobal,xi,yi)

Extract the bathymetry from the NetCDF file `bathname`. The parameter `isglobal`
 is true if the NetCDF file covers the whole globe and thus the last longitude
point can be considered to be right next to the first longitude point.
`xi` and `yi` are vectors defining the bounding box of the data. No
interpolation is performed.

**Convention:** b is positive in the water and negative in the air.
"""
function extract_bath(bath_name,isglobal,xi,yi)

    info("Extracting bathymetry from file: $(bath_name)")
    if isglobal == true
        info("Working with a global bathymetry");
    end;

    dxi = xi[2] - xi[1]
    dyi = yi[2] - yi[1]

    nc = Dataset(bath_name)
    x = nc["lon"].var[:]
    y = nc["lat"].var[:]

    rx = x[2]-x[1];
    ry = y[2]-y[1];
    X0 = x[1];
    Y0 = y[1];

    redx = dxi/rx;
    redy = dyi/ry;

    i0 = round(Int,(xi[1]-dxi-X0)/rx)+1;
    i1 = round(Int,(xi[end]+dxi-X0)/rx)+1;
    i=i0:i1;

    j0 = round(Int,(yi[1]-dyi-Y0)/ry)+1;
    j1 = round(Int,(yi[end]+dyi-Y0)/ry)+1;
    j=j0:j1;

    b = zeros(length(i),length(j));

    if isglobal
        i2 = mod.(i-1,size(nc["bat"] ,1))+1;
        jumps = [0; find(abs.(i2[2:end]-i2[1:end-1]) .> 1); length(i2)];

        for l=1:length(jumps)-1
            b[(jumps[l]+1):jumps[l+1],:] = nc["bat"].var[i2[jumps[l]+1]:i2[jumps[l+1]],j];
        end
    else
        i = maximum([minimum(i) 1]):minimum([maximum(i) length(x)]);
        j = maximum([minimum(j) 1]):minimum([maximum(j) length(y)]);
        b[:,:] = nc["bat"].var[i,j];
    end

    bx = X0 + rx*(i-1)
    by = Y0 + ry*(j-1)

    return bx,by,-b
end

"""
    xi,yi,bath = divand.load_bath(bath_name,isglobal,xi,yi)

Load the bathymetry from the NetCDF file `bathname`. The parameter `isglobal` is true if the NetCDF file covers
the whole globe and thus the last longitude point can be considered to be right next to the first longitude point.
`xi` and `yi` are vectors containing the longitude and latitude grid onto which the bathymetry should be
interpolated.

"""

function load_bath(bath_name,isglobal,xi,yi)

    bx,by,b = extract_bath(bath_name,isglobal,xi,yi)

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
    Xi,Yi = ndgrid(xi,yi);

    itp = interpolate((bx,by), b, Gridded(Linear()))
    bi = itp[xi,yi];

    return  xi,yi,bi
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

function load_mask(bath_name,isglobal,xi,yi,level::Number)

    info("Creating land-sea mask on level: $(level)")

    bx,by,b = extract_bath(bath_name,isglobal,xi,yi)

    # hack
    #b = -b
    #b[isnan.(b)] = 10
    # end hack
    #@show extrema(b)
    mask = b .> level;


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
    Xi,Yi = ndgrid(xi,yi);

    itp = interpolate((bx,by), Int.(mask),Gridded(Linear()))
    mif = itp[xi,yi];

    mi = mif .> 1/2;

    return  xi,yi,mi
end


function load_mask(bath_name,isglobal,x,y,levels::AbstractVector)
    data = [load_mask(bath_name,isglobal,x,y,level) for level in levels ];

    return data[1][1],data[1][2],cat(3,[d[3] for d in data]...)
end
