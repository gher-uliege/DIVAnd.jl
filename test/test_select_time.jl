using Base.Test
include("../src/select_time.jl")

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

TS = TimeSelectorYW(years,yearwindow,monthlists)
@test length(TS) == 8


centraltimes = ctimes(TS)

@test length(centraltimes) == length(TS)

@test Dates.Year(centraltimes[1]).value == 1993

sel = select(TS,1,obstime)

@test all(years[1] - yearwindow/2 .<= Dates.year.(obstime[sel]) .<= years[1] + yearwindow/2)
@test all(1 .<= Dates.month.(obstime[sel]) .<= 3)



# running average
times = DateTime(1990,1,1):Dates.Month(1):DateTime(2010,12,31)
window = 90
TS = TimeSelectorRunningAverage(times,window)
@test length(TS) == length(times)
sel = select(TS,1,obstime)
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

TS = TimeSelectorYearListMonthList(yearlists,monthlists)
@test length(TS) == 8

centraltimes = ctimes(TS)

@test length(centraltimes) == length(TS)

@test Dates.Year(centraltimes[1]).value == 1993

sel = select(TS,1,obstime)

@test all(years[1] - yearwindow/2 .<= Dates.year.(obstime[sel]) .<= years[1] + yearwindow/2)
@test all(1 .<= Dates.month.(obstime[sel]) .<= 3)




monthlists = [
    [2,3],
    [4,5,6],
    [7,8,9],
    [10:12; 1]
];

TS = TimeSelectorYearListMonthList(yearlists,monthlists)
@test length(TS) == 8

centraltimes = ctimes(TS)

@test length(centraltimes) == length(TS)

@test Dates.Year(centraltimes[1]).value == 1993

sel = select(TS,4,obstime)

@test all(years[1] - yearwindow/2 .<= Dates.year.(obstime[sel]) .<= years[1] + yearwindow/2)
obsmonth = Dates.month.(obstime[sel])

@test all( (10 .<= obsmonth .<= 12) .| (obsmonth .== 1))
