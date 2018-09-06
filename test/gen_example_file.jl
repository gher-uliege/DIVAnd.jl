using NCDatasets
if VERSION >= v"0.7"
    using Random
    using Dates
end

nobs = 1000;

randrg(min,max,nobs) = min .+ (max-min) * rand(Float32,nobs)

obslon = randrg(3.0418334f0, 11.81f0, nobs)
obslat = randrg(42.0f0, 44.0f0, nobs)
obsdepth = randrg(-0.0f0, 2762.0f0, nobs)
obstime = rand(DateTime(1892,09,25):Dates.Day(1):DateTime(2017,10,02),nobs)

Salinity = Float32.(30 .+ obsdepth / 3000 + 2 .* obslon/10 + 2 .* obslat/20)

obsid = repeat(Vector{Char}("486-1451364"),inner = (1,nobs))

filename = joinpath(dirname(@__FILE__),"..","data","sample-file.nc")

ds = Dataset(filename,"c")
# Dimensions

ds.dim["observations"] = nobs
ds.dim["idlen"] = size(obsid,1)

# Declare variables

ncobslon = defVar(ds,"obslon", Float32, ("observations",))
ncobslon.attrib["units"] = "degrees_east"
ncobslon.attrib["standard_name"] = "longitude"
ncobslon.attrib["long_name"] = "longitude"

ncobslat = defVar(ds,"obslat", Float32, ("observations",))
ncobslat.attrib["units"] = "degrees_north"
ncobslat.attrib["standard_name"] = "latitude"
ncobslat.attrib["long_name"] = "latitude"

ncobstime = defVar(ds,"obstime", Float64, ("observations",))
ncobstime.attrib["units"] = "days since 1900-01-01 00:00:00"
ncobstime.attrib["standard_name"] = "time"
ncobstime.attrib["long_name"] = "time"

ncobsdepth = defVar(ds,"obsdepth", Float32, ("observations",))
ncobsdepth.attrib["units"] = "meters"
ncobsdepth.attrib["positive"] = "down"
ncobsdepth.attrib["standard_name"] = "depth"
ncobsdepth.attrib["long_name"] = "depth below sea level"

ncobsid = defVar(ds,"obsid", Char, ("idlen", "observations"))
ncobsid.attrib["long_name"] = "observation identifier"
ncobsid.attrib["coordinates"] = "obstime obsdepth obslat obslon"

ncSalinity = defVar(ds,"Salinity", Float32, ("observations",))

ncobslon[:] = obslon
ncobslat[:] = obslat
ncobstime[:] = obstime
ncobsdepth[:] = obsdepth
ncobsid[:] = obsid
ncSalinity[:] = Salinity

close(ds)
