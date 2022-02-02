using Test
using DelimitedFiles
using Statistics
using Random
using StableRNGs
import DIVAnd

rng = Random.GLOBAL_RNG
rng = StableRNG(1234)

# test data for basic statistics
x = [1.0, 2.0, 3.0, 4.0]
y = -[1.0, 2.0, 3.0, 4.0]
sumx = sum(x)
sumx2 = sum(x .^ 2)

sumy = sum(y)
sumy2 = sum(y .^ 2)

sumxy = sum(x .* y)

meanx, stdx = DIVAnd.stats(sumx, sumx2, length(x))
meanx, meany, stdx, stdy, covar, corr =
    DIVAnd.stats(sumx, sumx2, sumy, sumy2, sumxy, length(x))

@test meanx ≈ mean(x)
@test stdx ≈ std(x)

@test corr ≈ -1

# DIVA fit

mincount = 50

mask, (pm, pn), (xi, yi) = DIVAnd.DIVAnd_squaredom(2, range(0, stop = 1, length = 100))

lenx = 0.05;
leny = 0.05;
Nens = 1
distbin = 0:0.02:0.3
mincount = 100

Random.seed!(rng,1234)
field = DIVAnd.random(mask, (pm, pn), (lenx, leny), Nens, rng = rng)

x = (xi[:], yi[:])
v = field[:]

distx, covar, corr, varx, count, stdcovar = DIVAnd.empiriccovarmean(x, v, distbin, mincount)

minlen = 0.001
maxlen = 0.1

var0opt = covar[1]
L = range(minlen, stop = maxlen, length = 100);

mu, K, len_scale = DIVAnd.DIVAnd_kernel(2, [1, 2, 1])

J(L) = sum(((covar - var0opt * K.(distx * len_scale / L)) ./ stdcovar) .^ 2)
Jmin, imin = findmin(J.(L))
lenopt = L[imin]


var0opt, lensopt, distx, covar, fitcovar = DIVAnd.fit_isotropic(x, v, distbin, mincount)


@test lensopt ≈ lenx rtol = 0.2
@test var0opt ≈ 1 rtol = 0.5

# port of DIVA fit from Fortran

#fname = "/home/abarth/src/DIVA/DIVA3D/src/Fortran/Util/smalltestdata.txt"
#A = readdlm(fname) :: Array{Float64,2}

#const fitlen = DIVAnd.fitlen

A = readdlm(joinpath(dirname(@__FILE__), "..", "data", "testdata.txt"))::Array{Float64,2}

n = size(A, 1)
x = A[:, 1]
y = A[:, 2]
d = A[:, 3]
weight = A[:, 4]

# use all pairs
nsamp = 0
varbak, RL, dbinfo = DIVAnd.fitlen((x, y), d, weight, nsamp; rng = rng)

# reference value are  from DIVA fit (Fortran version)
# git commit
# cb243004ffca6b49797f53dc3ccc357d71759cd1 (Fri Mar 30 16:48:10 2018 +0200)

@test 1.43710339 ≈ RL rtol = 1e-6
@test 1.39825165 ≈ dbinfo[:sn] rtol = 1e-6
@test 24.1159477 ≈ varbak rtol = 1e-6
@test 0.75906783342361450 ≈ dbinfo[:rqual] rtol = 1e-6

# random samples
nsamp = 150
varbak, RL, dbinfo = DIVAnd.fitlen((x, y), d, weight, nsamp; rng = rng)

# reference value from Julia implementation with seed set to 150
# fluctuations are large for different seeds

@test 1.6184793426171522 ≈ RL
@test 1.3058876287893113 ≈ dbinfo[:sn]
@test 23.425097759577064 ≈ varbak
@test 0.6885893124253655 ≈ dbinfo[:rqual]

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


mask, (pm, pn, po), (xi, yi, zi) =
    DIVAnd.DIVAnd_squaredom(3, range(0, stop = 1, length = 30))

lenx = leny = lenz = 0.2
Nens = 1
Random.seed!(rng,12345)
field = DIVAnd.random(mask, (pm, pn, po), (lenx, leny, lenz), Nens, rng = rng)


z = [0.3, 0.5, 0.7]
s = 1:7:length(field)
x = (xi[s], yi[s], zi[s])
v = field[s]
epsilon2 = ones(length(x[3])) + x[3][:] .^ 2

fitlenxy,
dbinfo = @test_logs (:info, r".*at*") match_mode = :any DIVAnd.fithorzlen(
    x,
    v,
    z;
    epsilon2 = epsilon2,
    rng = rng,
);

@test median(fitlenxy) ≈ lenx rtol = 0.3

s = 1:3:length(field)
x = (xi[s], yi[s], zi[s])
v = field[s]
epsilon2 = ones(length(x[3])) + x[3][:] .^ 2

fitlenz,
dbinfo = @test_logs (:info, r".*at*") match_mode = :any DIVAnd.fitvertlen(
    x,
    v,
    z;
    epsilon2 = epsilon2,
    rng = rng,
);
@test median(fitlenz) ≈ lenz rtol = 0.5


# RandomCoupels iterators

iseed = 1234
Random.seed!(rng,iseed)

@test collect(DIVAnd.RandomCoupels(1000, 10, rng)) ==
      collect(DIVAnd.RandomCoupels(1000, 10, rng))
