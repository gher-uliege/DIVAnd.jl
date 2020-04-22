using DIVAnd
using Dates
using Random

# Warm-up grid
dx, dy = 0.5, 0.5
lonr0 = 11.5:dx:20
latr0 = 39:dy:46
depthr0 = collect(0.:15.:30.);
monthlist0 = [[1]];

# Final grid
dx, dy = 0.125, 0.125
lonr = 11.5:dx:20
latr = 39:dy:46
timerange = [Date(1950,1,1),Date(2017,12,31)];
depthr = [0., 5., 10., 20., 30., 50., 100., 200., 300., 500., 750., 1000.];
varname = "Salinity"
yearlist = [1900:2017];
monthlist = [[1,2,3], [4,5,6], [7,8,9], [10,11,12]];

TS0 = DIVAnd.TimeSelectorYearListMonthList(yearlist, monthlist0);
TS = DIVAnd.TimeSelectorYearListMonthList(yearlist, monthlist);
@show TS;

datadir = "./Adriatic/"
datafile = joinpath(datadir, "AdriaticSea_SDC.txt")
if !isdir(datadir)
    @info("Creating data directory")
    mkdir(datadir)
end

if !isfile(datafile)
    @info("Downloading data file")
    download("https://dox.ulg.ac.be/index.php/s/A4Eu9nEoovYLtGr/download", datafile)
else
    @info("Data file already downloaded")
end

@time obsval,obslon,obslat,obsdepth,obstime,obsid = ODVspreadsheet.load(Float64,[datafile],
                           ["Water body salinity"]; nametype = :localname );

bathname = "./Adriatic/gebco_30sec_4.nc"
if !isfile(bathname)
    download("https://dox.ulg.ac.be/index.php/s/U0pqyXhcQrXjEUX/download",bathname)
else
    @info("Bathymetry file already downloaded")
end

@time bx0, by0, b0 = load_bath(bathname, true, lonr0, latr0);

mask0 = falses(size(b0,1),size(b0,2),length(depthr0))
for k = 1:length(depthr0)
    for j = 1:size(b0, 2)
        for i = 1:size(b0, 1)
            mask0[i,j,k] = b0[i,j] >= depthr0[k]
        end
    end
end

@time bx, by, b = load_bath(bathname,true,lonr,latr);

mask1 = falses(size(b,1),size(b,2),length(depthr))
for k = 1:length(depthr)
    for j = 1:size(b,2)
        for i = 1:size(b,1)
            mask1[i,j,k] = b[i,j] >= depthr[k]
        end
    end
end
@show size(mask1)


sz0 = (length(lonr0), length(latr0), length(depthr0));
lenx = fill(50_000.,sz0)   # 100 km
leny = fill(50_000.,sz0)   # 100 km
lenz = fill(0.,sz0);
len0 = (lenx, leny, lenz);
epsilon2 = 0.1;

sz = (length(lonr),length(latr),length(depthr));
lenx = fill(50_000.,sz)   # 100 km
leny = fill(50_000.,sz)   # 100 km
lenz = fill(0.,sz);
len1 = (lenx, leny, lenz);
epsilon2 = 0.1;

filename0 = "Water_body_$(replace(varname," "=>"_"))_$(randstring(['a':'z'; '0':'9'], 12)).nc"

if isfile(filename0)
    rm(filename0) # delete the previous analysis
    @info "Removing file $filename01"
end

@info("Warm-up run")
@time dbinfo = diva3d((lonr0, latr0, depthr0, TS0),
    (obslon,obslat,obsdepth,obstime), obsval,
    len0, epsilon2,
    filename0, varname,
    bathname=bathname,
    mask = mask0,
    fitcorrlen = false,
    niter_e = 2
    );

filename0 = "Water_body_$(replace(varname," "=>"_"))_$(randstring(['a':'z'; '0':'9'], 12)).nc"

if isfile(filename0)
    rm(filename0) # delete the previous analysis
    @info "Removing file $filename01"
end

@info("Race run")
@time dbinfo = diva3d((lonr, latr, depthr, TS),
(obslon,obslat,obsdepth,obstime), obsval,
len1, epsilon2,
filename0, varname,
bathname=bathname,
mask = mask1,
fitcorrlen = false,
niter_e = 2
);
