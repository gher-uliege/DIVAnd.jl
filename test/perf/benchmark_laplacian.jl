using Base.Cartesian
using Base.Test
import divand

sz = (1000,1000)
#sz = (100,100)
#sz = (10,10)

x = randn(sz)

ii = [Float64(i) for i = 1:sz[1], j = 1:sz[2]]
jj = [Float64(j) for i = 1:sz[1], j = 1:sz[2]]


x = [Float64(i^2) for i = 1:sz[1], j = 1:sz[2]]

pmn = (2*ones(sz),ones(sz))
nu = (ones(sz),ones(sz))

#pmn = (2*ones(sz),ones(sz))

#x = [Float64(i^2 + j^2/3) for i = 1:sz[1], j = 1:sz[2]]
pmn = (ii,jj+ii/10)
nu = (ii+2jj,jj+3ii)


function test_sp(mask,pmn,nu,x0,Nmax=1)
    sz = size(x0)
    x = x0[mask]

    L = divand.divand_laplacian(Val{:sparse},mask,pmn,nu,[false,false])

    @time @inbounds for nt = 1:Nmax
        x = L*x
    end

    Lx = zeros(x0)
    Lx[mask] = x
    return Lx
end


function test_sp_inplace(mask,pmn,nu,x0,Nmax=1)
    sz = size(x0)
    x = x0[mask]

    L = divand.divand_laplacian(Val{:sparse},mask,pmn,nu,[false,false])
    Lx = similar(x)

    @time @inbounds for nt = 1:Nmax
        #x = L*x
        A_mul_B!(Lx,L,x)
        x[:] = Lx
    end

    Lx2 = zeros(x0)
    Lx2[mask] = x
    return Lx2
end
    

function test_lap8(mask,pmn,nu,x0,Nmax=1)
    ivol,nus = divand.divand_laplacian_prepare(mask,pmn,nu)

    x = copy(x0)
    Lx = similar(x)
    
    @time @inbounds for nt = 1:Nmax
        divand.divand_laplacian_apply!(ivol,nus,x,Lx)
        x[:] = Lx
    end

    return x
end


mask = trues(sz)
mask[3:4,3:4] = false

Nmax = 10

Lxsp0 = test_sp(mask,pmn,nu,x,Nmax);
include("../../src/override_ssmult.jl")
Lxsp = test_sp(mask,pmn,nu,x,Nmax);
Lx2 = test_lap8(mask,pmn,nu,x,Nmax);
Lxsp2 = test_sp_inplace(mask,pmn,nu,x,Nmax);


@test Lxsp ≈ Lx2 atol=1e-4

