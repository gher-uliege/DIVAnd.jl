

struct AbstractTimeSelector
end




"""
    TS = TimeSelectorYearListMonthList(yearlists,monthlists)

The structure `TS` handels the time arregation base on `yearlists` and 
`monthlists`. `yearlists` is vector of ranges  (containing start and end year), 
for example `[1980:1990,1990:2000,2000:2010]`.

`monthlists` is vector of two-element vector (containing start and end month), for 
example `[1:3,4:6,7:9,10:12]`

If a month range spans beyond December, then all Months must be specified, e.g.
example `[2:4,5:6,7:9,[10,11,12,1]]` or `[2:4,5:6,7:9,[10:12;1]]`. 
However using `[2:4,5:6,7:9,10:1]` (bug!) will result in
an empty month range.
"""


struct TimeSelectorYearListMonthList
    yearlists
    monthlists
end

Base.length(TS::TimeSelectorYearListMonthList) = length(TS.yearlists) * length(TS.monthlists)

function ctimes(TS::TimeSelectorYearListMonthList)
    timeclim = DateTime[]

    for yearrange in TS.yearlists
        @assert(length(yearrange) > 0)
        yearc = yearrange[end ÷ 2]
        
        for monthrange in TS.monthlists
            @assert(length(monthrange) > 0)

            # central time instance
            timecentral =
                # day 16 of the central months
                if length(monthrange) % 2 == 1
                    DateTime(yearc,monthrange[(end+1) ÷ 2],16,0,0,0)
                else
                    DateTime(yearc,monthrange[end÷2 + 1],1,0,0,0)
                end
            
            push!(timeclim,timecentral)
        end
    end
    
    return timeclim
end

function select(TS::TimeSelectorYearListMonthList,index,obstime)   
    yearindex = (index-1) ÷ length(TS.monthlists) +1
    mlindex = (index-1) % length(TS.monthlists) +1

    yearlist = TS.yearlists[yearindex]
    monthlist = TS.monthlists[mlindex]

    # convertion to Int is necessary on 32-bit systems
    s = falses(Int.(size(obstime)))

    # loop over all observation time instance
    for i = 1:length(obstime)
        # s[i] is true if the observation is within the time range
        s[i] = yearlist[1] <= Dates.year(obstime[i]) <= yearlist[end]

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


struct TimeSelectorRunningAverage
    times # central times
    window # in days
end

Base.length(TS::TimeSelectorRunningAverage) = length(TS.times)
ctimes(TS::TimeSelectorRunningAverage) = TS.times

function select(TS::TimeSelectorRunningAverage,index,obstime)   
    # convertion to Int is necessary on 32-bit systems
    s = falses(Int.(size(obstime)))

    # loop over all observation time instance
    for i = 1:length(obstime)
        s[i] =  Dates.Millisecond(obstime[i] - TS.times[index]).value <= 1000*24*60*60*TS.window
    end

    return s
end

"""
    TS = TimeSelectorYW(years,yearwindow,monthlists)

The structure `TS` handels the time arregation base on `years` and 
`monthlists`. It is similar to `TimeSelectorYearListMonthList` except that 
the elements of `yearlists` are centred around `years` and span 
`yearwindow` years. `yearlists` is in fact constructed by additing and subtract 
`yearwindow/2` to every  element of years.

"""
function TimeSelectorYW(years,yearwindow,monthlists)
    yearlists = [y-yearwindow/2:y+yearwindow/2 for y in years]
    return TimeSelectorYearListMonthList(yearlists,monthlists)
end

