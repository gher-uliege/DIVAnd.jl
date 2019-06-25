if VERSION >= v"0.7.0-beta.0"
    using Test
    using Compat
    using Random
    using Statistics
    using Dates
else
    using Base.Test
    using Compat: range
end
using DIVAnd

#include("../src/utils.jl")

# gradient

x1,x2 = DIVAnd.ndgrid([1:2:20;],[1:3:30;])
pm = ones(size(x1))/2
pn = ones(size(x1))/3
f = 2. * x1 + x2

fx,fy = DIVAnd.cgradient((pm,pn),f)

@test fx ≈ fill(2,size(f))
@test fy ≈ fill(1,size(f))


x1,x2,x3 = DIVAnd.ndgrid([1:2:20;],[1:3:30;],[1:2:20;])
pm = ones(size(x1))/2
pn = ones(size(x1))/3
po = ones(size(x1))/2

f = 2. * x1 + x2 + 4. * x3

fx,fy,fz = DIVAnd.cgradientn((pm,pn,po),f)

@test fx ≈ fill(2,size(f))
@test fy ≈ fill(1,size(f))
@test fz ≈ fill(4,size(f))


# fill

c = randn(10,10)
valex = -9999
c[3:5,6:10] .= valex
cf = DIVAnd.ufill(c,valex);
@test sum(cf == valex) == 0

c = randn(10,10,20)
valex = -9999
c[3:5,6:10,1:4] .= valex
cf = DIVAnd.ufill(c,valex);
@test sum(cf == valex) == 0


# lengraddepth

mask,(pm,pn),(xi,yi) = DIVAnd.DIVAnd_squaredom(
    2,range(-10, stop = 10, length = 100))
h = 1000 * (tanh.(xi) .+ 1);
L = 2.
RL = DIVAnd.lengraddepth((pm,pn),h,L)

@test maximum(RL) < 1 + 10*eps(Float64)
@test RL[40,1] < RL[60,1]



# Greens functions for 1D diffusion
# 1/sqrt(4 π k t) * exp(-x^2 / (4kt))

z = range(-50,stop = 50,length = 201);
f = zeros(size(z))
f[(end+1)÷2] = 1

filterscale = 10
ff = DIVAnd.smoothfilter(z,f,filterscale)

fref = (z[2]-z[1]) * exp.(-z.^2/(2*filterscale^2)) / sqrt(2* π * filterscale^2)

@test sum(fref) ≈ 1 atol=1e-4
@test sum(ff) ≈ 1 atol=1e-4
@test maximum(abs.(ff - fref)) < 1e-4

#clf(); plot(z,ff, label = "sol"); plot(z,fref,label = "ref"); legend()

# random field
mask,(pm,pn),(xi,yi) = DIVAnd_rectdom(
    range(0, stop = 1, length = 100),
    range(0, stop = 1, length = 110))

lenx = .05;
leny = .05;
Nens = 100
pmn = (pm,pn)
len = (lenx,leny)
if VERSION >= v"0.7.0-beta.0"
   Random.seed!(123)
else
   srand(123)
end
field = DIVAnd.random(mask,pmn,len,Nens)
@test size(field) == (size(mask,1),size(mask,2),Nens)

@test std(field[50,50,:]) ≈ 1 atol=0.2

# interpolation

x1 = collect(1.:2:20)
x2 = collect(1.:3:30)
X1,X2 = DIVAnd.ndgrid(x1,x2)
f = 2. * X1 + X2

@test DIVAnd.interp((X1,X2),f,([3.],[1.])) ≈ [7.]
@test DIVAnd.interp((x1,x2),f,([3.],[1.])) ≈ [7.]


# weigts of observations

# two observation close-by should have less weight than the 3rd observation
# far away
weight = DIVAnd.weight_RtimesOne(([0.,0.1,2],[0.,0.,0.]),[1.,1.])
@test weight[1] < weight[3]


# days since

@test DIVAnd.dayssince(DateTime(1900,1,1); t0 = DateTime(1900,1,1)) == 0
@test DIVAnd.dayssince(DateTime(1900,1,2); t0 = DateTime(1900,1,1)) == 1



# flood-fill

mask = trues(8,8)
mask[3,:] .= false

m = DIVAnd.floodfillpoint(mask,CartesianIndex(1,1))
@test all(m[1:2,:])
@test all(.!m[3:end,:])

index = DIVAnd.floodfill(mask)
@test index[1,1] == 2
@test index[end,1] == 1
