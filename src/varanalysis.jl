"""
work1, work2: size of mask

"""

function Bsqrt!{T}(sv,coeff,ivol,nus,nmax,α,x::Array{T,1},work1,work2,Bsqrtx)
    work2[:] = 0
    work2[sv.mask[1]] = x
    work2[:] = work2[:] .* sqrt.(ivol[:])

    for niter = 1:(nmax ÷ 2)
        divand_laplacian_apply!(ivol,nus,work2,work1)
        #work2 += α * work1
        BLAS.axpy!(α,work1,work2)
    end

    work2[:] = coeff .* work2

    Bsqrtx[:] = pack(sv,(work2,))
end


"""
Variational analysis similar to 3D-var

Kernel is the solution of the n-dimensional diffusion equation

∂c/∂t =  ∇ ⋅ (D ∇ c)

n-dimensional Green’s function

G(x,x',t) = (4πDt)^(-n/2)  exp( - |x -x'|² / (4Dt))
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
    α0 = 1/(2 * n * max([maximum(pmn[i].^2 .* nu[i]) for i in 1:n]...)) :: T
    
    # 10% safety margin
    α = α0 / 1.1

    # number of iterations
    #nmax = round(Int,1/(2*α))
    # number of iterations 1/(2*α) (round to the closest be even number)
    nmax = 2*round(Int,1/(4*α))
    @show nmax

    # the background error covariance matrix is
    # B =  (4π α nmax)^(n/2) prod(Ld) (I + α * D)^nmax;

    # x^T B x is the integral which takes also the volumn of each grid cell into account

    # W is the norm
    # xᵀ W y correspond to the integral ∫ f g dx

    ivol,nus = divand.divand_laplacian_prepare(mask,pmn,nu)

    vol = 1.0 ./ ivol
    
    #W = Diagonal(vol)

    coeff = ((4π * α * nmax)^(n/4) * sqrt(prod(Ld))) * sqrt.(vol)
   
    H = H * Diagonal(sqrt.(ivol[mask]))
    
    work1 = zeros(size(mask))
    work2 = zeros(size(mask))
    tmpx = zeros(s.sv.n)    
    b = zeros(s.sv.n)
                        
    function fun!(x,fx)
        # tmpx = B^1/2 x
        Bsqrt!(s.sv,coeff,ivol,nus,nmax,α,x,work1,work2,tmpx)        
        # tmpx = B^1/2 * H' * (R \ (H * B^1/2 x ))
        Bsqrt!(s.sv,coeff,ivol,nus,nmax,α,H' * (R \ (H * tmpx)),work1,work2,tmpx)
        # fx = x + B^1/2 * H' * (R \ (H * B^1/2 x ))
        fx[:] = x + tmpx
    end

    # b = B^{1/2} * H' * (R \ yo)
    Bsqrt!(s.sv,coeff,ivol,nus,nmax,α,(H' * (R \ yo)),work1,work2,b)

    #@code_warntype Bsqrt(n,s.sv,coeff,ivol,nus,nmax,α,H' * (R \ yo),work1,work2)

    # adjust tolerance
    tol = tol * s.sv.n / length(yo)

    # xp = (I + B^1/2 * H' * (R^{-1} * (H * B^1/2)))^{-1} b

    #@show divand.checksym(s.sv.n,fun!)
          
    xp,success,s.niter = divand.conjugategradient(
        fun!,b; tol = tol, maxit = maxit,
        progress = progress);

    # tmpx = B^1/2 * xp
    Bsqrt!(s.sv,coeff,ivol,nus,nmax,α,xp,work1,work2,tmpx)
    
    return unpack(s.sv, Diagonal(sqrt.(ivol[mask])) * tmpx,NaN)[1],s
end
