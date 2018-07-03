using Base.Cartesian
import DIVAnd
if VERSION >= v"0.7.0-beta.0"
    using Test
else
    using Base.Test
end


function test_sp(mask,pmn,nu,x0,Nmax=1)
    sz = size(x0)
    x = x0[mask]

    L = DIVAnd.DIVAnd_laplacian(Val{:sparse},mask,pmn,nu,falses(ndims(mask)))

    @inbounds for nt = 1:Nmax
        x = L*x
    end

    Lx = zeros(size(x0))
    Lx[mask] = x
    return Lx
end


function test_lap8(mask,pmn,nu,x0,Nmax=1)
    ivol,nus = DIVAnd.DIVAnd_laplacian_prepare(mask,pmn,nu)
    x = copy(x0)

    @inbounds for nt = 1:Nmax
        x = DIVAnd.DIVAnd_laplacian_apply(ivol,nus,x)
    end

    return x
end

#sz = (1000,1000)  # too large for a 32-bit OS
Nmax = 1

for sz in [(20,),(100,100),(20,20,20),(5,5,5,5)]
    x = randn(sz)
    mask = trues(sz)

    ij = DIVAnd.ndgrid([Float64.(1:s) for s in sz]...)
    x = ij[1].^2

    if length(sz) == 2
        pmn = (ij[1],ij[2]+ij[1]/10)
        nu = (ij[1]+2 * ij[2],ij[2]+3 * ij[1])

        mask[3:4,3:4] .= false
    else
        pmn = ntuple(i -> ones(sz),length(sz))
        nu = ntuple(i -> i*ones(sz),length(sz))
    end



    Lxsp = test_sp(mask,pmn,nu,x,Nmax)
    Lx2 = test_lap8(mask,pmn,nu,x,Nmax)

    ivol,nus = DIVAnd.DIVAnd_laplacian_prepare(mask,pmn,nu)

    Lx = DIVAnd.DIVAnd_laplacian_apply(ivol,nus,x)

    @test Lxsp ≈ Lx2

    # check symmetry
    D = DIVAnd.DIVAnd_laplacian(Val{:sparse},mask,pmn,nu,falses(ndims(mask)));

    vol = 1 ./ .*(pmn...);
    D2 = DIVAnd.sparse_diag(vol[mask]) * D;
    @test maximum(abs.(D2 - D2')) < 1e-9
end

# degenerated case
sz = (10,1)
mask = trues(sz)
x = zeros(sz)
pmn = (zeros(sz),zeros(sz))
nu = (zeros(sz),zeros(sz))
for j = 1:sz[2]
    for i = 1:sz[1]
        nu[1][i,j] = i + 2*j
        nu[2][i,j] = j + 3*i

        pmn[1][i,j] = i
        pmn[2][i,j] = j+i/10

        x[i,j] = i^2
    end
end

Lxsp = test_sp(mask,pmn,nu,x,Nmax)
ivol,nus = DIVAnd.DIVAnd_laplacian_prepare(mask,pmn,nu)
Lx = DIVAnd.DIVAnd_laplacian_apply(ivol,nus,x)
@test Lxsp ≈ Lx
