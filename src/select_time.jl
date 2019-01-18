

struct AbstractTimeSelector
end




"""
    TS = TimeSelectorYearListMonthList(yearlists,monthlists)

The structure `TS` handles the time aggregation based on `yearlists` and
`monthlists`. `yearlists` is a vector of ranges (containing start and end years),
for example `[1980:1990,1990:2000,2000:2010]`.

`monthlists` is a vector of two-element vector (containing start and end months), for
example `[1:3,4:6,7:9,10:12]`

If a month range spans beyond December, then all Months must be specified, e.g.
example `[2:4,5:6,7:9,[10,11,12,1]]` or `[2:4,5:6,7:9,[10:12;1]]`.
However using `[2:4,5:6,7:9,10:1]` (bug!) will result in
an empty month range.

## Example

# seasonal climatology using all data from 1900 to 2017
# for winter (December-February), spring, summer, autumn

TS = DIVAnd.TimeSelectorYearListMonthList([1900:2017],[[12,1,2],[3,4,5],[6,7,8],[9,10,11]])

"""
struct TimeSelectorYearListMonthList{T1 <: AbstractVector,T2 <: AbstractVector}
    yearlists::T1
    monthlists::T2
end

Base.length(TS::TimeSelectorYearListMonthList) = length(TS.yearlists) * length(TS.monthlists)

function ctimes(TS::TimeSelectorYearListMonthList)
    timeclim = DateTime[]

    for yearrange in TS.yearlists
        @assert(length(yearrange) > 0)
        yearc = yearrange[(end+1) ÷ 2]

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


"""
    timecim = timestart(TS::TimeSelectorYearListMonthList)

Return the start date of all intervals defined by `TS`.
"""
function timesstart(TS::TimeSelectorYearListMonthList)
    timeclim = DateTime[]

    for yearrange in TS.yearlists
        @assert(length(yearrange) > 0)

        for monthrange in TS.monthlists
            @assert(length(monthrange) > 0)

            # start time instance
            time0 = DateTime(yearrange[1],monthrange[1],1,0,0,0)
            push!(timeclim,time0)
        end
    end

    return timeclim
end

"""
    timecim = timeend(TS::TimeSelectorYearListMonthList)

Return the end date of all intervals defined by `TS`.
"""
function timesend(TS::TimeSelectorYearListMonthList)
    timeclim = DateTime[]

    for yearrange in TS.yearlists
        @assert(length(yearrange) > 0)

        for monthrange in TS.monthlists
            @assert(length(monthrange) > 0)

            # end time instance
            time0 = Dates.lastdayofmonth(DateTime(yearrange[end],monthrange[end],1,0,0,0))
            push!(timeclim,time0)
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

"""
    TS = TimeSelectorRunningAverage(times,window)

The structure `TS` handles the time aggregation based on vector of
central `times` and the time window given in days.
Observations at the i-th time instance will be selected
if the dates is between `times[i]-w/2` and `time[i]+w/2` where
`w` is the time window expressed as days.
"""
struct TimeSelectorRunningAverage{T1 <: AbstractVector, T2 <: Number}
    times::T1 # central times
    window::T2 # in days
end

Base.length(TS::TimeSelectorRunningAverage) = length(TS.times)

ctimes(TS::TimeSelectorRunningAverage) = TS.times

timesstart(TS::TimeSelectorRunningAverage) =
    TS.times - Dates.Millisecond(round(Int64,Int64(TS.window) * 24*60*60*1000/2))

timesend(TS::TimeSelectorRunningAverage) =
    TS.times + Dates.Millisecond(round(Int64,Int64(TS.window) * 24*60*60*1000/2))

function select(TS::TimeSelectorRunningAverage,index,obstime)
    # convertion to Int is necessary on 32-bit systems
    s = falses(Int.(size(obstime)))

    # loop over all observation time instance
    for i = 1:length(obstime)
        s[i] =  abs(Dates.Millisecond(obstime[i] - TS.times[index]).value) <= 1000*24*60*60*TS.window
    end

    return s
end

"""
    TS = TimeSelectorYW(years,yearwindow,monthlists)

The structure `TS` handles the time aggregation based on `years` and
`monthlists`. It is similar to `TimeSelectorYearListMonthList` except that
the elements of `yearlists` are centred around `years` and span
`yearwindow` years. `yearlists` is in fact constructed by adding and subtracting
`yearwindow/2` to every element of years.

"""
function TimeSelectorYW(years,yearwindow,monthlists)
    yearlists = [y-yearwindow/2:y+yearwindow/2 for y in years]
    return TimeSelectorYearListMonthList(yearlists,monthlists)
end



"""
    cbounds = climatology_bounds(TS)

Produce an matrix for DateTimes where `cbounds[i,1]` is the start time of
all sub-intervals defined by this `i`-th time instance and
`cbounds[i,2]` is the end time of all sub-intervals defined by this `i`-th time
instance.
"""
function climatology_bounds(TS)
# https://web.archive.org/web/20180326074452/http://cfconventions.org/Data/cf-conventions/cf-conventions-1.7/cf-conventions.html

# it has a climatology attribute, which names a variable with
# dimensions (n,2), n being the dimension of the climatological time
# axis. Using the units and calendar of the time coordinate variable,
# element (i,0) of the climatology variable specifies the beginning of
# the first subinterval and element (i,1) the end of the last
# subinterval used to evaluate the climatological statistics with
# index i in the time dimension.

    b = Array{DateTime}(undef,2,length(TS))
    b[1,:] = timesstart(TS)
    b[2,:] = timesend(TS)
    return b
end
