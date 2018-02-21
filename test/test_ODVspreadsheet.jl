using Base.Test
import divand
ODV = divand.ODVspreadsheet


fname = joinpath(dirname(@__FILE__),"..","data","sample_ODV.txt")

odv = ODV.readODVspreadsheet(fname)

T = Float64
data,data_qv,lon,lat,depth,depth_qv,time,time_qv,EDMO,LOCAL_CDI_ID = ODV.loadprofile(T,odv,1,"SDN:P01::PSSTTS01")

@test EDMO[1] == "1234"

data,data_qv,lon,lat,depth,depth_qv,time,time_qv,EDMO,LOCAL_CDI_ID = ODV.loadprofile(T,odv,1,"SDN:P01::SLCAAAZX")

@test data[1] == 37.5
@test EDMO[1] == "1234"

fnames = [joinpath(dirname(@__FILE__),"..","data",n) for n in ["sample_ODV.txt","sample_ODV2.txt"]]

P01names = ["SDN:P01::ODSDM021"]
profiles,lons,lats,depths,times,ids = ODV.load(T,fnames,P01names)

@test 30 in profiles
@test 31 in profiles

fname = joinpath(dirname(@__FILE__),"..","data","sample_ODV_aggregated.txt")
datanames = ["Water body salinity"]

profiles,lons,lats,depths,times,ids = ODV.load(T,[fname],datanames,nametype = :localname)

@test length(profiles) > 0

# test if data with bad depth information are discarded

fname_qv = joinpath(dirname(@__FILE__),"..","data","sample_ODV_qv.txt")
profiles,lons,lats,depths,times,ids = ODV.load(T,[fname_qv],["SDN:P01::SLCAAAZX"])

@test 4. in depths # quality flag "good - 1"
@test !(3. in depths) # quality flag "missing value - 9"

# values from
# https://web.archive.org/web/20171129142108/https://www.hermetic.ch/cal_stud/chron_jdate.htm
# rounded to 3 hour

@test divand.ODVspreadsheet.parsejd(2454142.125) == DateTime(2007,02,10,03,0,0)

# values from
# http://www.julian-date.com/ (setting GMT offset to zero)
# https://web.archive.org/web/20180212213256/http://www.julian-date.com/

@test divand.ODVspreadsheet.parsejd(2455512.375) == DateTime(2010,11,11,9,0,0)
