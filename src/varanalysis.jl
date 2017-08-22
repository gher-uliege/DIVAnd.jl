function Bsqrt{T}(n,sv,ivol,nus,Ld,nmax,α,x::Array{T,1},Lx)
    #@code_warntype unpack(sv,x)
    
    xup = unpack(sv,x)[1]
    
    for niter = 1:(nmax ÷ 2)
        divand_laplacian_apply!(ivol,nus,xup,Lx)
        #xup += α * Lx
        BLAS.axpy!(α,Lx,xup)
    end

    xup = (4π * α * nmax)^(n/4) * sqrt(prod(Ld)) * (sqrt.(ivol) .* xup)

    #return pack(sv,(xup,))[:,1]::Array{T,1}
    return pack(sv,(xup,))
end

len_harmonize{T <: Number,N}(len::T,mask::AbstractArray{Bool,N})::NTuple{N, Array{T,N}} = ((fill(len,size(mask)) for i=1:N)...)
len_harmonize{T <: Number,N}(len::NTuple{N,T},mask::AbstractArray{Bool,N})::NTuple{N, Array{T,N}} = ((fill(len[i],size(mask)) for i=1:N)...)
function len_harmonize{T <: Number,N}(len::NTuple{N,AbstractArray{T,N}},mask::AbstractArray{Bool,N})::NTuple{N, Array{T,N}}
    # for i=1:N
    #     if size(mask) != size(len[i])
    #         error("mask ($(formatsize(size(mask)))) and correlation length ($(formatsize(size(len[i])))) have incompatible size")
    #     end
    # end

    return len
end




"""
Variational analysis similar to 3D-var

Kernel is the solution of the n-dimensional diffusion equation

∂c/∂t =  ∇ ⋅ (D ∇ c)

n-dimensional Green’s function

G(x,x',t) = (4πDt)^(-n/2)  exp( - |x -x'|² / (4Dt))
http://www.rpgroup.caltech.edu/~natsirt/aph162/diffusion_old.pdf

"""

function varanalysis{T,N}(mask::AbstractArray{Bool,N},pmn,xi,x,f::AbstractVector{T},len,epsilon2; tol::T = 1e-5)
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

    dx_min = 1/max([maximum(pm) for pm in pmn]...) :: T
    nu_max = max([maximum(ν) for ν in nu]...) :: T
    # need to proof this (currently this is just an analogy based on 2D)
    α0 = dx_min^2/(2 * nu_max * n)

    # 10% safety margin
    α = α0 / 1.1

    # number of iterations
    nmax = round(Int,1/(2*α))

    # the background error covariance matrix is
    # B =  (4π α nmax)^(n/2) prod(Ld) (I + α * D)^nmax;

    # x^T B x is the integral which takes also the volumn of each grid cell into account 

    ivol,nus = divand.divand_laplacian_prepare(mask,pmn,nu)

    Lx = zeros(size(mask))
    
    function fun!(x,fx)
        fx[:] = x + Bsqrt(n,s.sv,ivol,nus,Ld,nmax,α,H' * (R \ (H * (Bsqrt(n,s.sv,ivol,nus,Ld,nmax,α,x,Lx)))),Lx)
    end

    b = Bsqrt(n,s.sv,ivol,nus,Ld,nmax,α,(H' * (R \ yo)),Lx)

    #@code_warntype Bsqrt(n,s.sv,ivol,nus,Ld,nmax,α,H' * (R \ yo),Lx)

    # adjust tolerance
    tol = tol * s.sv.n / length(yo)

    xp,success,s.niter = divand.conjugategradient(fun!,b; tol = tol);
    #@code_warntype divand.conjugategradient(fun!,b; tol = tol);

    xa = Bsqrt(n,s.sv,ivol,nus,Ld,nmax,α,xp,Lx)
    
    return unpack(s.sv,xa,NaN)[1],s
end
