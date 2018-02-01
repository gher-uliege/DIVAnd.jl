using Base.Test
import divand
#include("/home/abarth/projects/Julia/divand.jl/src/utils.jl")

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
