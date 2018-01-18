using Base.Test


struct TimeSelectorRunningAverage
    times # central times
    window # in days
end

Base.length(TS::TimeSelectorRunningAverage) = length(TS.times)
ctimes(TS::TimeSelectorRunningAverage) = TS.times

function select(TS::TimeSelectorRunningAverage,index,obstime)   
    s = falses(size(obstime))

    # loop over all observation time instance
    for i = 1:length(obstime)
        s[i] =  Dates.Millisecond(obstime[i] - TS.times[index]).value <= 1000*24*60*60*TS.window
    end

    return s
end




struct TimeSelectorYS
    years
    monthlists
    year_window
end

Base.length(TS::TimeSelectorYS) = length(TS.years)*length(TS.monthlists)

function ctimes(TS::TimeSelectorYS)
    timeclim = DateTime[]

    for yearc in TS.years
        for monthlist in TS.monthlists
            # central time instance
            timecentral =
                # day 16 of the central months
                if length(monthlist) % 2 == 1
                    DateTime(yearc,monthlist[(end+1) ÷ 2],16,0,0,0)
                else
                    DateTime(yearc,monthlist[end÷2 + 1],1,0,0,0)
                end
            
            push!(timeclim,timecentral)
        end
    end
    
    return timeclim
end

function select(TS::TimeSelectorYS,index,obstime)
    yearindex = (index-1) ÷ length(TS.monthlists) +1
    mlindex = (index-1) % length(TS.monthlists) +1

    yearc = TS.years[yearindex]
    yearrange = [yearc - year_window/2, yearc + year_window/2]
    monthlist = TS.monthlists[mlindex]
    
    s = falses(size(obstime))

    @show yearc,monthlist
    # loop over all observation time instance
    for i = 1:length(obstime)
        # s[i] is true if the observation is within the time range
        s[i] = yearrange[1] <= Dates.year(obstime[i]) <= yearrange[end]

        # sm is true if the month is one of the months is monthlist
        sm = false
        for m in monthlist
            sm = sm || (Dates.month(obstime[i]) == m)
        end

        # keep an observation is year and month are suitable
        s[i] = s[i] && sm;
    end
    return s    
end






years = 1993:1994 
year_window = 10

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

TS = TimeSelectorYS(years,monthlists,year_window)
@test length(TS) == 8


centraltimes = ctimes(TS)

@test length(centraltimes) == length(TS)

@test Dates.Year(centraltimes[1]).value == 1993

sel = select(TS,1,obstime)

@test all(years[1] - year_window/2 .<= Dates.year.(obstime[sel]) .<= years[1] + year_window/2)
@test all(1 .<= Dates.month.(obstime[sel]) .<= 3)



obstime = DateTime(1990,1,1) : DateTime(2010,12,31)

# running average
times = DateTime(1990,1,1):Dates.Month(1):DateTime(2010,12,31)
window = 90
TS = TimeSelectorRunningAverage(times,window)
@test length(TS) == length(times)
sel = select(TS,1,obstime)
@test all((obstime[sel] - times[1]) .<= Dates.Day(window))





#------------

"""

`yearlists` is vector of two-element vector (containing start and end year), for 
example [[1980,1990],[1990,2000],[2000,2010]]


`monthlists` is vector of two-element vector (containing start and end month), for 
example [[2,7],[8,1]]

"""
struct TimeSelectorYearListMonthList
    yearlists # list of 
    monthslist # in days
end

Base.length(TS::TimeSelectorYearListMonthList) = length(TS.yearlist) * length(TS.monthlist)

function ctimes(TS::TimeSelectorYearListMonthList)
    timeclim = DateTime[]

    for yearrange in TS.yearlists
        for monthrange in TS.monthlists

            if monthrange[2] < monthrange[1]
                # add 12 month
                monthrange[2] += 12
            end

            @assert(monthrange[2] >= monthrange[1])

            # range 
            ml = monthrange[1]:monthrange[2]
                
            # central time instance
            timecentral =
                # day 16 of the central months
                if length(ml) % 2 == 1
                    # we must make sure that the central month is between
                    # 1 and 12
                    centralmonth = (ml[(end+1) ÷ 2] -1) % 12 +1
                    DateTime(yearc,centralmonth,16,0,0,0)
                else
                    centralmonth = (ml[end÷2 + 1] -1) % 12 +1
                    DateTime(yearc,centralmonth,1,0,0,0)
                end
            
            push!(timeclim,timecentral)
        end
    end
    
    return timeclim
end

function select(TS::TimeSelectorYearListMonthList,index,obstime)   
    s = falses(size(obstime))

    # loop over all observation time instance
    for i = 1:length(obstime)
        s[i] =  Dates.Millisecond(obstime[i] - TS.times[index]).value <= 1000*24*60*60*TS.window
    end

    return s
end

