using Base.Test
import divand

#include("../src/utils.jl")

# gradient

x1,x2 = divand.ndgrid([1:2:20;],[1:3:30;])
pm = ones(size(x1))/2
pn = ones(size(x1))/3
f = 2. * x1 + x2

fx,fy = divand.cgradient((pm,pn),f)

@test fx ≈ fill(2,size(f))
@test fy ≈ fill(1,size(f))


x1,x2,x3 = divand.ndgrid([1:2:20;],[1:3:30;],[1:2:20;])
pm = ones(size(x1))/2
pn = ones(size(x1))/3
po = ones(size(x1))/2

f = 2. * x1 + x2 + 4. * x3

fx,fy,fz = divand.cgradientn((pm,pn,po),f)

@test fx ≈ fill(2,size(f))
@test fy ≈ fill(1,size(f))
@test fz ≈ fill(4,size(f))


# fill

c = randn(10,10)
valex = -9999
c[3:5,6:10] = valex
cf = divand.ufill(c,valex);
@test sum(cf == valex) == 0

c = randn(10,10,20)
valex = -9999
c[3:5,6:10,1:4] = valex
cf = divand.ufill(c,valex);
@test sum(cf == valex) == 0


# lengraddepth

mask,(pm,pn),(xi,yi) = divand.divand_squaredom(2,linspace(-10,10,100))
h = 1000 * (tanh.(xi)+1);
L = 2.
RL = divand.lengraddepth((pm,pn),h,L)

@test maximum(RL) < 1 + 10*eps(Float64)
@test RL[40,1] < RL[60,1]



# Greens functions for 1D diffusion
# 1/sqrt(4 π k t) * exp(-x^2 / (4kt))

z = linspace(-50,50,201);
f = zeros(z)
f[(end+1)÷2] = 1

scale = 10
ff = divand.smoothfilter(z,f,scale)

fref =  (z[2]-z[1]) * exp.(-z.^2/(2*scale^2)) / sqrt(2* π * scale^2);

@test sum(fref) ≈ 1 atol=1e-4
@test sum(ff) ≈ 1 atol=1e-4
@test maximum(abs.(ff - fref)) < 1e-4

#clf(); plot(z,ff, label = "sol"); plot(z,fref,label = "ref"); legend()


# random field


xi,yi = divand.ndgrid(linspace(0,1,100),linspace(0,1,110))

mask = trues(size(xi))
pm = ones(size(xi)) / (xi[2,1]-xi[1,1])
pn = ones(size(xi)) / (yi[1,2]-yi[1,1])
lenx = .05;
leny = .05;
Nens = 100
pmn = (pm,pn)
len = (lenx,leny)
srand(123)
field = divand.random(mask,pmn,len,Nens)
@test size(field) == (size(mask,1),size(mask,2),Nens)

@test std(field[50,50,:]) ≈ 1 atol=0.2

# interpolation

x1 = collect(1.:2:20)
x2 = collect(1.:3:30)
X1,X2 = divand.ndgrid(x1,x2)
f = 2. * X1 + X2

@test divand.interp((X1,X2),f,([3.],[1.])) ≈ [7.]
@test divand.interp((x1,x2),f,([3.],[1.])) ≈ [7.]


# weigts of observations

# two observation close-by should have less weight than the 3rd observation
# far away
weight = divand.weight_RtimesOne(([0.,0.1,2],[0.,0.,0.]),[1.,1.])
@test weight[1] < weight[3]


# days since

@test divand.dayssince(DateTime(1900,1,1); t0 = DateTime(1900,1,1)) == 0
@test divand.dayssince(DateTime(1900,1,2); t0 = DateTime(1900,1,1)) == 1
