using NCDatasets
using Missings
using DIVAnd
using Compat: @info, range

NCSDN = DIVAnd.NCSDN


basedir = joinpath(dirname(@__FILE__),"..","..",
                   "DIVAnd-example-data","SDN_D85_WP8_netCDF_files_examples")

fname = joinpath(basedir,"netCDF_vertical_profiles_ctd.nc")

param = "TEMPPR01"

T = Float64
fnames = [fname]

data,lon,lat,z,time,ids = NCSDN.load(T,fnames,param;
     qualityflags = [NCSDN.GOOD_VALUE, NCSDN.PROBABLY_GOOD_VALUE])

fnames = [
    joinpath(basedir,"netCDF_timeseries_tidegauge.nc"),
    joinpath(basedir,"netCDF_timeseries_tidegauge_with_instrument.nc"),
    joinpath(basedir,"netCDF_trajectory_meteorological_data.nc"),
    joinpath(basedir,"netCDF_trajectory_tsg_with_instrument.nc"),
    joinpath(basedir,"netCDF_vertical_profiles_ctd.nc"),
    joinpath(basedir,"netCDF_vertical_profiles_ctd_with_instruments.nc"),
    joinpath(basedir,"netCDF_vertical_profiles_xbt_with_fall_rate_and_instruments.nc")]

data,lon,lat,z,time,ids = NCSDN.load(T,fnames,param;
     qualityflags = [NCSDN.GOOD_VALUE, NCSDN.PROBABLY_GOOD_VALUE])


#param = "TEMPET01"

#data,lon,lat,z,time,ids = NCSDN.load(T,fnames,param;
#     qualityflags = [NCSDN.GOOD_VALUE, NCSDN.PROBABLY_GOOD_VALUE])


nothing