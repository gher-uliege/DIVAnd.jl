# download data to run test script

basedir = joinpath(dirname(@__FILE__),"..","..","divand-example-data")
bathname = joinpath(basedir,
                    "Global","Bathymetry","gebco_30sec_16.nc")

obsname = joinpath(basedir,
                   "Provencal","WOD-Salinity.nc")


mkpath(joinpath(basedir,"Global","Bathymetry"))
mkpath(joinpath(basedir,"Provencal"))

if !isfile(bathname)
    download("https://b2drop.eudat.eu/s/o0vinoQutAC7eb0/download",bathname)
end

if !isfile(obsname)
    download("https://b2drop.eudat.eu/s/UsF3RyU3xB1UM2o/download",obsname)
end

