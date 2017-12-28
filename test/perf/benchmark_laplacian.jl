using Base.Cartesian
using Base.Test
import divand


function test_sp(mask,pmn,nu,x0,Nmax=1)
    sz = size(x0)
    x = x0[mask]

    L = divand.divand_laplacian(Val{:sparse},mask,pmn,nu,falses(ndims(mask)))

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

    L = divand.divand_laplacian(Val{:sparse},mask,pmn,nu,falses(ndims(mask)))
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

#sz = (1000,1000)  # too large for a 32-bit OS
Nmax = 10

#for sz in [(20,),(100,100),(20,20,20),(5,5,5,5)]
sz = (200,200,5)

    x = randn(sz)
    mask = trues(sz)

    ij = divand.ndgrid([Float64.(1:s) for s in sz]...)
    x = ij[1].^2

    if length(sz) == 2
        pmn = (ij[1],ij[2]+ij[1]/10)
        nu = (ij[1]+2 * ij[2],ij[2]+3 * ij[1])

        mask[3:4,3:4] = false        
    else
        pmn = ntuple(i -> ones(sz),length(sz))
        nu = ntuple(i -> i*ones(sz),length(sz))
    end




Lxsp0 = test_sp(mask,pmn,nu,x,Nmax);
include("../../src/override_ssmult.jl")
Lxsp = test_sp(mask,pmn,nu,x,Nmax);
Lx2 = test_lap8(mask,pmn,nu,x,Nmax);
Lxsp2 = test_sp_inplace(mask,pmn,nu,x,Nmax);


@test Lxsp ≈ Lx2 atol=1e-4

nothing
