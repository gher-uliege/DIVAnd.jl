using Base.Test
import divand
ODV = divand.ODVspreadsheet


fname = joinpath(dirname(@__FILE__),"..","data","sample_ODV.txt")

odv = ODV.readODVspreadsheet(fname)

T = Float64
data,data_qv,lon,lat,depth,time,time_qv,EDMO,LOCAL_CDI_ID = ODV.loadprofile(T,odv,1,"SDN:P01::PSSTTS01")

@test EDMO[1] == "1234"

fnames = [joinpath(dirname(@__FILE__),"..","data",n) for n in ["sample_ODV.txt","sample_ODV2.txt"]]
@show fnames

P01names = ["SDN:P01::ODSDM021"]
profiles,lons,lats,depths,times,ids = ODV.load(T,fnames,P01names)

@test 30 in profiles
@test 31 in profiles

fname = joinpath(dirname(@__FILE__),"..","data","sample_ODV_aggregated.txt")
datanames = ["Water body salinity"]

profiles,lons,lats,depths,times,ids = ODV.load(T,[fname],datanames,nametype = :localname)

@test length(profiles) > 0
