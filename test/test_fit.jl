using Base.Test
import divand

# test data for basic statistics
x = [1.,2.,3.,4.]
y = -[1.,2.,3.,4.]
sumx = sum(x)
sumx2 = sum(x.^2)

sumy = sum(y)
sumy2 = sum(y.^2)

sumxy = sum(x .* y)

meanx,stdx = divand.stats(sumx,sumx2,length(x))
meanx,meany,stdx,stdy,covar,corr = divand.stats(sumx,sumx2,sumy,sumy2,sumxy,length(x))

@test meanx ≈ mean(x)
@test stdx ≈ std(x)

@test corr ≈ -1


# DIVA fit

mincount = 50

xi,yi = divand.ndgrid(linspace(0,1,100),linspace(0,1,100))
mask = trues(size(xi))
pm = ones(size(xi)) / (xi[2,1]-xi[1,1])
pn = ones(size(xi)) / (yi[1,2]-yi[1,1])

lenx = .05;
leny = .05;
Nens = 1
distbin = 0:0.05:1
mincount = 100

srand(1234)
field = divand.random(mask,(pm,pn),(lenx,leny),Nens)

x = (xi[:],yi[:])
v = field[:]

distx,covar,corr,varx,count = divand.empiriccovar(
    x,v,distbin,mincount)


var0opt,lensopt,distx,covar,fitcovar = divand.fit_isotropic(
    x,v,distbin,mincount)


@test lensopt ≈ lenx rtol=0.2
@test var0opt ≈ 1 rtol=0.5
