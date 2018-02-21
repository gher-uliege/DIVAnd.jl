using Base.Test
import divand

obstime = DateTime(1990,1,1) : DateTime(2010,12,31)




years = 1993:1994 
yearwindow = 10

# winter: January-March    1,2,3
# spring: April-June       4,5,6
# summer: July-September   7,8,9
# autumn: October-December 10,11,12

monthlists = [
    [1,2,3],
    [4,5,6],
    [7,8,9],
    [10,11,12]
];

TS = divand.TimeSelectorYW(years,yearwindow,monthlists)
@test length(TS) == 8


centraltimes = divand.ctimes(TS)

@test length(centraltimes) == length(TS)

@test Dates.Year(centraltimes[1]).value == 1993

sel = divand.select(TS,1,obstime)

@test all(years[1] - yearwindow/2 .<= Dates.year.(obstime[sel]) .<= years[1] + yearwindow/2)
@test all(1 .<= Dates.month.(obstime[sel]) .<= 3)



# running average
times = DateTime(1990,1,1):Dates.Month(1):DateTime(2010,12,31)
window = 90
TS = divand.TimeSelectorRunningAverage(times,window)
@test length(TS) == length(times)
sel = divand.select(TS,1,obstime)
@test all((obstime[sel] - times[1]) .<= Dates.Day(window))








#-------------------------

yearwindow = 10
yearlists = [y-yearwindow/2:y+yearwindow/2 for y in 1993:1994]

# winter: January-March    1,2,3
# spring: April-June       4,5,6
# summer: July-September   7,8,9
# autumn: October-December 10,11,12

monthlists = [
    [1,2,3],
    [4,5,6],
    [7,8,9],
    [10,11,12]
];

TS = divand.TimeSelectorYearListMonthList(yearlists,monthlists)
@test length(TS) == 8

centraltimes = divand.ctimes(TS)

@test length(centraltimes) == length(TS)

@test Dates.Year(centraltimes[1]).value == 1993

sel = divand.select(TS,1,obstime)

@test all(years[1] - yearwindow/2 .<= Dates.year.(obstime[sel]) .<= years[1] + yearwindow/2)
@test all(1 .<= Dates.month.(obstime[sel]) .<= 3)




monthlists = [
    [2,3],
    [4,5,6],
    [7,8,9],
    [10:12; 1]
];

TS = divand.TimeSelectorYearListMonthList(yearlists,monthlists)
@test length(TS) == 8

centraltimes = divand.ctimes(TS)

@test length(centraltimes) == divand.length(TS)

@test Dates.Year(centraltimes[1]).value == 1993

sel = divand.select(TS,4,obstime)

@test all(years[1] - yearwindow/2 .<= Dates.year.(obstime[sel]) .<= years[1] + yearwindow/2)
obsmonth = Dates.month.(obstime[sel])

@test all( (10 .<= obsmonth .<= 12) .| (obsmonth .== 1))
