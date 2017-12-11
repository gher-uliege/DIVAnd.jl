using Base.Test

#import divand

#ODVspreadsheet = divand.ODVspreadsheet

import ODVspreadsheet

fname = joinpath(dirname(@__FILE__),"..","data","sample_ODV.txt")


odv = ODVspreadsheet.readODVspreadsheet(fname)

T = Float64
data,data_qv,lon,lat,depth,time,time_qv,EDMO,LOCAL_CDI_ID = ODVspreadsheet.loadprofile(T,odv,1,"SDN:P01::PSSTTS01")

@test EDMO[1] == "1234"


data,data_qv,lon,lat,depth,time,time_qv,EDMO,LOCAL_CDI_ID = ODVspreadsheet.loadprofile(T,odv,1,"SDN:P01::PSSTTS01")

fnames = [joinpath(dirname(@__FILE__),"..","data",n) for n in ["sample_ODV.txt","sample_ODV2.txt"]]

P01names = ["SDN:P01::ODSDM021"]
profiles,lons,lats,depths,times,ids = ODVspreadsheet.load(T,fnames,P01names)
