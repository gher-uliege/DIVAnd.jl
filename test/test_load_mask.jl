using Test
using DIVAnd

bathname = joinpath(
    dirname(@__FILE__),
    "..",
    "..",
    "DIVAnd-example-data",
    "Global",
    "Bathymetry",
    "gebco_30sec_16.nc",
)

bathname8 = joinpath(
    dirname(@__FILE__),
    "..",
    "..",
    "DIVAnd-example-data",
    "Global",
    "Bathymetry",
    "gebco_30sec_8.nc",
)

if !isfile(bathname)
    @info("download bathymetry $bathname")
    bathname = download("https://dox.ulg.ac.be/index.php/s/U0pqyXhcQrXjEUX/download")
end

if !isfile(bathname8)
    @info("download bathymetry $bathname8")
    bathname8 = download("https://dox.ulg.ac.be/index.php/s/wS6Y8P8NhIF60eG/download")
end

bathisglobal = true;


x0 = 27.0
x1 = 41.8
dx = 0.12
y0 = 40.3
y1 = 46.8
dy = 0.1
level = 0

lonr = x0:dx:x1
latr = y0:dy:y1

xi, yi, mi = load_mask(bathname, bathisglobal, lonr, latr, level)


levels = [0, 100, 1000]
xi, yi, mis = load_mask(bathname, bathisglobal, lonr, latr, levels)

@test mi == mis[:, :, 1]


lonr = -80:1:80
latr = -9:1:9
bx,by,b = load_bath(bathname,bathisglobal,lonr,latr);

# issue 45
dx = dy = 0.005
lonr = 19.45:dx:20.675
latr = 59.5:dy:59.9;

bx,by,b = load_bath(bathname8,bathisglobal,lonr,latr)
@test size(b) == (length(lonr),length(latr));


# global with wrapping of longitude
lonr = -180:1:180
latr = -89:1:89
bx,by,b = load_bath(bathname,bathisglobal,lonr,latr);
@test size(b) == (length(lonr),length(latr));

