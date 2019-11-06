using DIVAnd
using Statistics
using DelimitedFiles
using LinearAlgebra
using Random
using Test

Random.seed!(1)
NX1D = 10
LX1D = 8
xleft1D = -4
dx1D = LX1D / (NX1D)

xg1D = xleft1D + dx1D / 2:dx1D:xleft1D + LX1D
mask1D, pm1D, xi1D = DIVAnd.DIVAnd_rectdom(xg1D)



xo = randn(15)
inflation = ones(size(xo))

dens1D, LHM, LCV, LSCV = DIVAnd_heatmap(
    mask1D,
    pm1D,
    xi1D,
    (xo,),
    inflation,
    0;
    Ladaptiveiterations = 0,
)
#@show dens1D[2]
@test dens1D[2] ≈ 0.03390607052933598
dens1D, LHM, LCV, LSCV = DIVAnd_heatmap(
    mask1D,
    pm1D,
    xi1D,
    (xo,),
    inflation,
    1;
    Ladaptiveiterations = 0,
)
#@show dens1D[2]
@test dens1D[2] ≈ 0.053735329940045436
dens1D1, LHM, LCV, LSCV = DIVAnd_heatmap(
    mask1D,
    pm1D,
    xi1D,
    (xo,),
    inflation,
    1;
    Ladaptiveiterations = 1,
)
#@show dens1D1[2]
@test dens1D1[2] ≈ 0.03222559084580536
dens1D1, LHM, LCV, LSCV = DIVAnd_heatmap(
    mask1D,
    pm1D,
    xi1D,
    (xo,),
    inflation,
    1;
    Ladaptiveiterations = 2,
)
#@show dens1D1[2]
@test dens1D1[2] ≈ 0.038056634649593586
dens1D1nopt, LHM, LCV, LSCV = DIVAnd_heatmap(
    mask1D,
    pm1D,
    xi1D,
    (xo,),
    inflation,
    1;
    Ladaptiveiterations = 2,
    optimizeheat = false,
)
#@show dens1D1nopt[2]
@test dens1D1nopt[2] ≈ 0.03805663463970936
dens1D, LHM, LCV, LSCV = DIVAnd_heatmap(
    mask1D,
    pm1D,
    xi1D,
    (xo,),
    inflation,
    1;
    Ladaptiveiterations = 0,
    myheatmapmethod = "DataKernel",
)
#@show dens1D[2]
@test dens1D[2] ≈ 0.053735329940045436
dens1D, LHM, LCV, LSCV = DIVAnd_heatmap(
    mask1D,
    pm1D,
    xi1D,
    (xo,),
    inflation,
    1;
    Ladaptiveiterations = 0,
    myheatmapmethod = "GridKernel",
)
#@show dens1D[2]
@test dens1D[2] ≈ 0.05661676688382721
dens1D, LHM, LCV, LSCV = DIVAnd_heatmap(
    mask1D,
    pm1D,
    xi1D,
    (xo,),
    inflation,
    1;
    Ladaptiveiterations = 0,
    myheatmapmethod = "Automatic",
)
#@show dens1D[2]
@test dens1D[2] ≈ 0.05661676688382721


Random.seed!(1)
xo = randn(150)
inflation = ones(size(xo))

dens1D, LHM, LCV, LSCV = DIVAnd_heatmap(
    mask1D,
    pm1D,
    xi1D,
    (xo,),
    inflation,
    0;
    Ladaptiveiterations = 0,
    nmax = 100,
)
#@show dens1D[2]
@test dens1D[2] ≈ 0.033593831115271915


newcoord, newval, sumw, varp, idx = DIVAnd_superobs((xo,), inflation, 100)
@test newval[7] ≈ 1.0

newcoord, newval, sumw, varp, idx = DIVAnd_superobs(
    (xo,),
    inflation,
    100;
    intensive = false,
)
@test newval[7] ≈ 2.0
