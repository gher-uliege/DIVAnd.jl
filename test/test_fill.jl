using Base.Test


#imax = 100
#jmax = 100

c = randn(10,10)
valex = -9999
c[3:5,6:10] = valex
cf = ufill(c,valex);
@test sum(cf == valex) == 0

c = randn(10,10,20)
valex = -9999
c[3:5,6:10,1:4] = valex
cf = ufill(c,valex);
@test sum(cf == valex) == 0


#SSH = Dataset("/home/abarth/Utils/sossheig.nc")["sossheig"][:];
#c = copy(SSH.data); c[SSH.na] = valex;
#@time cf = ufill(c,valex);
#@show extrema(cf)

