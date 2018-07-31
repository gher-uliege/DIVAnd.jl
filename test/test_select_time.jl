if VERSION >= v"0.7.0-beta.0"
    using Test
else
    using Base.Test
end
import DIVAnd

obstime = DateTime(1990,1,1) : Dates.Day(1) : DateTime(2010,12,31)




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

TS = DIVAnd.TimeSelectorYW(years,yearwindow,monthlists)
@test length(TS) == 8


centraltime = DIVAnd.ctimes(TS)
@test Dates.Year(centraltime[1]).value == 1993

starttime = DIVAnd.timesstart(TS)
@test length(starttime) == length(TS)
@test Dates.Year(starttime[1]).value == 1993-yearwindow/2

sel = DIVAnd.select(TS,1,obstime)

@test all(years[1] - yearwindow/2 .<= Dates.year.(obstime[sel]) .<= years[1] + yearwindow/2)
@test all(1 .<= Dates.month.(obstime[sel]) .<= 3)



# running average
times = DateTime(1990,1,1):Dates.Month(1):DateTime(2010,12,31)
window = 90
TS = DIVAnd.TimeSelectorRunningAverage(times,window)
@test length(TS) == length(times)
sel = DIVAnd.select(TS,1,obstime)
@test all((obstime[sel] - times[1]) .<= Dates.Day(window))



centraltime = DIVAnd.ctimes(TS)
@test centraltime[1] == times[1]

starttime = DIVAnd.timesstart(TS)
@test starttime[1] == times[1]-Dates.Day(window/2)

endtime = DIVAnd.timesend(TS)
@test endtime[1] == times[1]+Dates.Day(window/2)





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

TS = DIVAnd.TimeSelectorYearListMonthList(yearlists,monthlists)
@test length(TS) == 8

starttime = DIVAnd.timesstart(TS)
endtime = DIVAnd.timesend(TS)

@test length(starttime) == length(TS)
@test length(endtime) == length(TS)

@test Dates.Year(starttime[1]).value == 1993-5
@test Dates.Year(endtime[1]).value == 1993+5

@test starttime[1] == DateTime(1993-5,1,1)
@test endtime[1] == DateTime(1993+5,3,31)

sel = DIVAnd.select(TS,1,obstime)

@test all(years[1] - yearwindow/2 .<= Dates.year.(obstime[sel]) .<= years[1] + yearwindow/2)
@test all(1 .<= Dates.month.(obstime[sel]) .<= 3)




monthlists = [
    [2,3],
    [4,5,6],
    [7,8,9],
    [10:12; 1]
];

TS = DIVAnd.TimeSelectorYearListMonthList(yearlists,monthlists)
@test length(TS) == 8

starttime = DIVAnd.timesstart(TS)

@test length(starttime) == DIVAnd.length(TS)

@test Dates.Year(starttime[1]).value == 1993-5

sel = DIVAnd.select(TS,4,obstime)

@test all(years[1] - yearwindow/2 .<= Dates.year.(obstime[sel]) .<= years[1] + yearwindow/2)
obsmonth = Dates.month.(obstime[sel])

@test all( (10 .<= obsmonth .<= 12) .| (obsmonth .== 1))
