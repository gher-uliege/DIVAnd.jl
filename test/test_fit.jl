using Base.Test
import divand
#include("../src/fit.jl")

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
distbin = 0:0.02:0.3
mincount = 100

srand(1234)
field = divand.random(mask,(pm,pn),(lenx,leny),Nens)


x = (xi[:],yi[:])
v = field[:]

distx,covar,corr,varx,count,stdcovar = divand.empiriccovarmean(
    x,v,distbin,mincount)

minlen = 0.001
maxlen = 0.1

var0opt = covar[1]
L = linspace(minlen,maxlen,100);

mu,K,len_scale = divand.divand_kernel(2,[1,2,1])

J(L) = sum(((covar - var0opt * K.(distx * len_scale/L))./stdcovar).^2)
Jmin,imin = findmin(J.(L))
lenopt = L[imin]


var0opt,lensopt,distx,covar,fitcovar = divand.fit_isotropic(
    x,v,distbin,mincount)


@test lensopt ≈ lenx rtol=0.2
@test var0opt ≈ 1 rtol=0.5

# port of DIVA fit from Fortran

#fname = "/home/abarth/src/DIVA/DIVA3D/src/Fortran/Util/smalltestdata.txt"
#A = readdlm(fname) :: Array{Float64,2}

#const fitlen = divand.fitlen

A = readdlm(joinpath(dirname(@__FILE__),"..","data","testdata.txt")) :: Array{Float64,2}

n = size(A,1)
x = A[:,1]
y = A[:,2]
d = A[:,3]
weight = A[:,4]

# use all pairs
nsamp = 0
varbak,RL,dbinfo = divand.fitlen((x,y),d,weight,nsamp)

# reference value are  from DIVA fit (Fortran version)
# git commit
# cb243004ffca6b49797f53dc3ccc357d71759cd1 (Fri Mar 30 16:48:10 2018 +0200)

@test 1.43710339 ≈ RL                       rtol=1e-6
@test 1.39825165 ≈ dbinfo[:sn]              rtol=1e-6
@test 24.1159477 ≈ varbak                   rtol=1e-6
@test 0.75906783342361450 ≈ dbinfo[:rqual]  rtol=1e-6

# random samples
nsamp = 150
varbak,RL,dbinfo = divand.fitlen((x,y),d,weight,nsamp)

# reference value from Julia implementation with seed set to 150
# fluctuations are large for different seeds

@test 2.6803941824646085  ≈ RL             rtol=0.1
@test 0.8076378487032164  ≈ dbinfo[:sn]    rtol=0.1
@test 18.48072398826336   ≈ varbak         rtol=0.1
@test 0.7338792203547216  ≈ dbinfo[:rqual] rtol=0.1

  
  
  
#=       
fname = "/home/abarth/src/DIVA/DIVA3D/src/Fortran/Util/testdata.txt"
A = readdlm(fname) :: Array{Float64,2}

n = size(A,1)
x = A[:,1]
y = A[:,2]
d = A[:,3]
weight = A[:,4]

# use all pairs
nsamp = 0
@time varbak,RL,dbinfo = fitlen((x,y),d,weight,nsamp)
@test 1.56002688           ≈ RL                       rtol=1e-6
@test 1.36453283           ≈ dbinfo[:sn]              rtol=1e-6
@test 25.4311714           ≈ varbak                   rtol=1e-6
@test 0.81123489141464233  ≈ dbinfo[:rqual]  rtol=1e-6
=#
