using NCDatasets
using Random
using Dates
using DataStructures

nobs = 1000;

randrg(min, max, nobs) = min .+ (max - min) * rand(Float32, nobs)

obslon = randrg(3.0418334f0, 11.81f0, nobs)
obslat = randrg(42.0f0, 44.0f0, nobs)
obsdepth = randrg(-0.0f0, 2762.0f0, nobs)
obstime = rand(DateTime(1892, 09, 25):Dates.Day(1):DateTime(2017, 10, 02), nobs)

Salinity = Float32.(30 .+ obsdepth / 3000 + 2 .* obslon / 10 + 2 .* obslat / 20)

obsid = repeat(Vector{Char}("486-1451364"), inner = (1, nobs))

filename = joinpath(dirname(@__FILE__), "..", "data", "sample-file.nc")

ds = Dataset(filename, "c")

# Dimensions

ds.dim["observations"] = nobs
ds.dim["idlen"] = size(obsid, 1)


# Declare variables

ncobslon = defVar(ds,"obslon", Float32, ("observations",), attrib = OrderedDict(
    "units"                     => "degrees_east",
    "standard_name"             => "longitude",
    "long_name"                 => "longitude",
))

ncobslat = defVar(ds,"obslat", Float32, ("observations",), attrib = OrderedDict(
    "units"                     => "degrees_north",
    "standard_name"             => "latitude",
    "long_name"                 => "latitude",
))

ncobstime = defVar(ds,"obstime", Float64, ("observations",), attrib = OrderedDict(
    "units"                     => "days since 1900-01-01 00:00:00",
    "standard_name"             => "time",
    "long_name"                 => "time",
))

ncobsdepth = defVar(ds,"obsdepth", Float32, ("observations",), attrib = OrderedDict(
    "units"                     => "meters",
    "positive"                  => "down",
    "standard_name"             => "depth",
    "long_name"                 => "depth below sea level",
))

ncobsid = defVar(ds,"obsid", Char, ("idlen", "observations"), attrib = OrderedDict(
    "long_name"                 => "observation identifier",
    "coordinates"               => "obstime obsdepth obslat obslon",
))

ncSalinity = defVar(ds,"Salinity", Float32, ("observations",))

ncobslon[:] = obslon
ncobslat[:] = obslat
ncobstime[:] = obstime
ncobsdepth[:] = obsdepth
ncobsid[:] = obsid
ncSalinity[:] = Salinity

close(ds)
