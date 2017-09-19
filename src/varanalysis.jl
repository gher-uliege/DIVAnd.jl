
function diffusion!(ivol,nus,α,nmax,x0,x)
    work1 = similar(x)
    x[:] = x0
    
    for niter = 1:nmax
        divand_laplacian_apply!(ivol,nus,x,work1)
        x[:] = x + α * work1        
    end

end
    

"""
work1, work2: size of mask

Symmetric matrix

SB = √(β) (1 + α L)^(nmax/2) W^{-1}

"""



function decompB!{T}(sv,β,ivol,nus,nmax,α,x::Array{T,1},work1,work2,decompBx)
    work2[:] = 0
    work2[sv.mask[1]] = x
    work2[:] = work2[:] .* ivol[:]

    for niter = 1:(nmax ÷ 2)
        divand_laplacian_apply!(ivol,nus,work2,work1)
        #work2 += α * work1
        BLAS.axpy!(α,work1,work2)
    end

    work2[:] = β * work2

    decompBx[:] = pack(sv,(work2,))
end


"""
Variational analysis similar to 3D-var

Kernel is the solution of the n-dimensional diffusion equation

∂c/∂t =  ∇ ⋅ (D ∇ c)

n-dimensional Green’s function

G(x,x',t) = (4πDt)^(-n/2)  exp( - |x -x'|² / (4Dt))


G(x,x',t) = det(D)^(-1/2) (4π t)^(-n/2)  exp( - (x -x')ᵀ D⁻¹ (x -x')ᵀ / (4t))

http://www.rpgroup.caltech.edu/~natsirt/aph162/diffusion_old.pdf

"""

function varanalysis{T,N}(mask::AbstractArray{Bool,N},pmn,xi,x,f::AbstractVector{T},len,epsilon2;
                          tol::T = 1e-5,
                          maxit::Int = 100000,
                          progress = (iter,x,r,tol2,fun!,b) -> nothing)
    
    n = ndims(mask)

    len = len_harmonize(len,mask)

    R = divand.divand_obscovar(epsilon2,length(f));

    s = divand.divand_struct(mask)

    # observation constrain
    constrain = divand.divand_obs(s,xi,x,f,R)
    yo = constrain.yo
    H = constrain.H

    Ld = T[mean(L) for L in len]
    nu = ([L.^2 for L in len]...) 

    # Building the Laplacian ∇ ⋅ (ν ∇ ϕ) where ν is the
    # correlation length-scale squared

    #dx_min = 1/max([maximum(pm) for pm in pmn]...) :: T
    #nu_max = max([maximum(ν) for ν in nu]...) :: T
    # need to proof this (currently this is just an analogy based on 2D)
    #α0 = dx_min^2/(2 * nu_max * n)
    #α0 = 1/(2 * n * max([maximum(pmn[i].^2 .* nu[i]) for i in 1:n]...)) :: T
    # ok
    α0 = 1/(2 * sum([maximum(pmn[i].^2 .* nu[i]) for i in 1:n])) :: T
    
    # 10% safety margin
    α = α0 / 1.1

    # number of iterations
    #nmax = round(Int,1/(2*α))
    # number of iterations 1/(2*α) (round to the closest be even number)
    nmax = 2*round(Int,1/(4*α))
    @show nmax

    # Green's functions
    # ∂c/∂t =  ∇ ⋅ (D ∇ c)

    # n-dimensional Green’s function
    # G(x,x',t) = (4πDt)^(-n/2)  exp( - |x -x'|² / (4Dt))
    # G(x,x',t) = det(D)^(-1/2) (4π t)^(-n/2)  exp( - (x -x')ᵀ D⁻¹ (x -x') / (4t))
    # 
    # c(x,t) = ∫ G(x,x',t) c₀(x') dx
    # 
    # In discrete where W is the norm (Δx_1 * Δx_2 * ... Δx_n)
    # (1 + α L)^nmax  x =  det(D)^(-1/2) (4π t)^(-n/2)  B W x
    
    
    # the background error covariance matrix is
    # B = β (I + α * D)^nmax W^{-1}
    # β = (4π α nmax)^(n/2) prod(Ld)
    # SB =  √(β) (I + α * D)^nmax/2 W^{-1}
    
    # B = SB W SB

    # W is the norm
    # xᵀ W y correspond to the integral ∫ f g dx

    ivol,nus = divand.divand_laplacian_prepare(mask,pmn,nu)

    vol = 1.0 ./ ivol

    sW = Diagonal(sqrt.(vol[mask]))
    
    #W = Diagonal(vol)

    @show Ld
    β = ((4π * α * nmax)^(n/4) * sqrt(prod(Ld)))   
    
    #H = H * Diagonal(sqrt.(ivol[mask]))
    
    work1 = zeros(size(mask))
    work2 = zeros(size(mask))
    tmpx = zeros(s.sv.n)    
    b = zeros(s.sv.n)

    # x + W^1/2 * SB^1/2 * H' * (R \ (H * SB^1/2 * W^1/2 * x ))
    function fun!(x,fx)
        # tmpx = SB^1/2 W^1/2 x
        decompB!(s.sv,β,ivol,nus,nmax,α,sW * x,work1,work2,tmpx)
        # tmpx = SB^1/2 * H' * (R \ (H * SB^1/2 x ))
        decompB!(s.sv,β,ivol,nus,nmax,α,H' * (R \ (H * tmpx)),work1,work2,tmpx)
        # fx = x + W^1/2 SB^1/2 * H' * (R \ (H * SB^1/2 x ))
        fx[:] = x + sW * tmpx
    end

    # b = W^1/2 SB^{1/2} * H' * (R \ yo)
    decompB!(s.sv,β,ivol,nus,nmax,α,(H' * (R \ yo)),work1,work2,b)
    b = sW * b 

    #@code_warntype decompB(n,s.sv,β,ivol,nus,nmax,α,H' * (R \ yo),work1,work2)

    # adjust tolerance
    tol = tol * s.sv.n / length(yo)

    # xp = (I + W^1/2 SB^1/2 * H' * (R^{-1} * (H * SB^1/2 * W^1/2)) )^{-1} b

    @show divand.checksym(s.sv.n,fun!)
          
    xp,success,s.niter = divand.conjugategradient(
        fun!,b; tol = tol, maxit = maxit,
        progress = progress);

    # tmpx = SB^1/2 * W^1/2 * xp
    decompB!(s.sv,β,ivol,nus,nmax,α,sW * xp,work1,work2,tmpx)

    @show "1",mean(β),nmax,α
    
    # debugging function
    function decompB(x)

        @show mean(β),nmax,α
        
        decompBx = similar(x)
        divand.decompB!(s.sv,β,ivol,nus,nmax,α,x,work1,work2,decompBx)

#        decompBx2 = (sW)^2 *  decompBx
#        divand.decompB!(s.sv,β,ivol,nus,nmax,α,decompBx2,work1,work2,decompBx)
        
        return decompBx        
    end

    info = Dict(:B => x -> decompB((sW)^2 * decompB(x)),:H => H, :R => R, :yo => yo,
                :mask => mask, :ivol => ivol, :nus => nus, :β => β,
                :α => α, :nmax => nmax)
    return unpack(s.sv, tmpx,NaN)[1],s,info
end
