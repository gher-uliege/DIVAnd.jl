using DIVAnd
using Test
using Dates

@test DIVAnd.unstagger([1, 2, 3]) == [0.5, 1.5, 2.5, 3.5]


@test DIVAnd.unstagger(1:3) == 0.5:1.0:3.5


x = 1:0.1:5

@test DIVAnd.findin(x, 2.1234) == 12
@test DIVAnd.findin(collect(x), 2.1234) == 12


mean_v,count = DIVAnd.binning((-1:1,-1:1),([0.5],[0.5]),[2.])

@test mean_v[end,end] == 2
@test count[end,end] == 1

@test all(in(Set([0,1])),count)



TS = DIVAnd.TimeSelectorYearListMonthList([1900:2020],[[12,1,2],[3,4,5],[6,7,8],[9,10,11]])

mean_v,count = DIVAnd.binning(
    (-1:1,-1:1,-1:1,TS),
    ([-1],[-1],[-1],[DateTime(1900,1,15)]),
    [2.])

@test mean_v[1,1,1,1] == 2
@test count[1,1,1,1] == 1
@test all(in(Set([0,1])),count)
