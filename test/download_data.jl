# download data to run test script

basedir = joinpath(dirname(@__FILE__),"..","..","DIVAnd-example-data")
bathname = joinpath(basedir,
                    "Global","Bathymetry","gebco_30sec_16.nc")

obsname = joinpath(basedir,
                   "Provencal","WOD-Salinity.nc")


mkpath(joinpath(basedir,"Global","Bathymetry"))
mkpath(joinpath(basedir,"Provencal"))

if !isfile(bathname)
    download("https://dox.ulg.ac.be/index.php/s/U0pqyXhcQrXjEUX/download",bathname)
end

if !isfile(obsname)
    download("https://dox.ulg.ac.be/index.php/s/PztJfSEnc8Cr3XN/download",obsname)
end

