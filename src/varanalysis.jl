
function diffusion!(ivol,nus,α,nmax,x0,x)
    work1 = similar(x)
    x[:] = x0

    for niter = 1:nmax
        DIVAnd_laplacian_apply!(ivol,nus,x,work1)
        x[:] = x + α * work1
    end

end


"""
work1, work2: size of mask

Symmetric matrix

SB = √(β) (1 + α L)^(nmax / 2) W^{-1}

where W is the volumne of the corresponding grid cell.
The background error covariance matrix B is SB W SB
"""
function decompB!(sv,β,ivol,nus,nmax,α,x::Array{T,1},work1,work2,decompBx) where T
    # does not work is some points are masked
    # for i in eachindex(work2)
    #     if sv.mask[1][i]
    #         work2[i] = x[i] * ivol[i]
    #     else
    #         work2[i] = 0
    #     end
    # end

    work2[:] .= 0
    work2[sv.mask[1]] = x
    work2[:] = work2[:] .* ivol[:]

    for niter = 1:(nmax ÷ 2)
        DIVAnd_laplacian_apply!(ivol,nus,work2,work1)
        #work2 += α * work1
        BLAS.axpy!(α,work1,work2)
    end

    for i in eachindex(work2)
        work2[i] = β * work2[i]
    end

    decompBx[:] = pack(sv,(work2,))
end

#function A_ldiv_B!(


"""
Variational analysis similar to 3D-var

Input:

  x0: start vector for iteration, at output it is the last state of the
   iteration. Note that x0 is related to the analysis xa by
      xa = SB^½ * W^½ * xa


  | x + W^½ * SB^½ * H' * (R \\ (H * SB^½ * W^½ * x ))   -   W^½ SB^{½} * H' * (R \\ yo) |
     <
  tol * s.sv.n / length(yo)  * | W^½ SB^{½} * H' * (R \\ yo) |

Kernel is the solution of the n-dimensional diffusion equation

∂c/∂t =  ∇ ⋅ (D ∇ c)

n-dimensional Green’s function

G(x,x',t) = (4πDt)^(-n/2)  exp( - |x -x'|² / (4Dt))


G(x,x',t) = det(D)^(-½) (4π t)^(-n/2)  exp( - (x -x')ᵀ D⁻¹ (x -x')ᵀ / (4t))

http://www.rpgroup.caltech.edu/~natsirt/aph162/diffusion_old.pdf

"""
function varanalysis(mask::AbstractArray{Bool,N},pmn,xi,x,
                     f::AbstractVector{T},len,epsilon2;
                     tol::T = 1e-5,
                     maxit::Int = 100000,
                     progress = (iter,x,r,tol2,fun!,b) -> nothing,
                     kwargs...) where N where T

    # all keyword agruments as a dictionary
    kw = Dict(kwargs)

    n = ndims(mask)

    len = len_harmonize(len,mask)

    R = DIVAnd.DIVAnd_obscovar(epsilon2,length(f));

    s = DIVAnd.DIVAnd_struct(mask)

    # observation constrain
    constrain = DIVAnd.DIVAnd_obs(s,xi,x,f,R)
    yo = constrain.yo
    H = constrain.H

    Ld = T[mean(L) for L in len]
    nu = ([L.^2 for L in len]...,)

    # Building the Laplacian ∇ ⋅ (ν ∇ ϕ) where ν is the
    # correlation length-scale squared

    #dx_min = 1/max([maximum(pm) for pm in pmn]...) :: T
    #nu_max = max([maximum(ν) for ν in nu]...) :: T
    # need to proof this (currently this is just an analogy based on 2D)
    #α0 = dx_min^2/(2 * nu_max * n)
    #α0 = 1/(2 * n * max([maximum(pmn[i].^2 .* nu[i]) for i in 1:n]...)) :: T
    # ok

    #@show [maximum(pmn[i].^2 .* nu[i]) for i in 1:n]
    α0 = 1/(2 * sum([maximum(pmn[i].^2 .* nu[i]) for i in 1:n])) :: T

    # 10% safety margin
    α = α0 / 1.1

    # number of iterations 1/(2*α) (round to the closest be even number)
    nmax = 2*round(Int,1/(4*α))
    #@show nmax

    # Green's functions
    # ∂c/∂t =  ∇ ⋅ (D ∇ c)

    # n-dimensional Green’s function
    # G(x,x',t) = (4πDt)^(-n/2)  exp( - |x -x'|² / (4Dt))
    # G(x,x',t) = det(D)^(-½) (4π t)^(-n/2)  exp( - (x -x')ᵀ D⁻¹ (x -x') / (4t))
    #
    # c(x,t) = ∫ G(x,x',t) c₀(x') dx
    #
    # In discrete where W is the norm (Δx_1 * Δx_2 * ... Δx_n)
    # (1 + α L)^nmax  x =  det(D)^(-½) (4π t)^(-n/2)  B W x


    # the background error covariance matrix is
    # B = β (I + α * D)^nmax W^{-1}
    # β = (4π α nmax)^(n/2) prod(Ld)
    # SB =  √(β) (I + α * D)^nmax/2 W^{-1}

    # B = SB W SB

    # W is the norm
    # xᵀ W y correspond to the integral ∫ f g dx

    ivol,nus = DIVAnd.DIVAnd_laplacian_prepare(mask,pmn,nu)

    vol = 1.0 ./ ivol

    sW = Diagonal(sqrt.(vol[mask]))

    #W = Diagonal(vol)

    #@show Ld
    β = ((4π * α * nmax)^(n/4) * sqrt(prod(Ld)))

    #H = H * Diagonal(sqrt.(ivol[mask]))

    work1 = zeros(size(mask))
    work2 = zeros(size(mask))
    tmpx = zeros(s.sv.n)
    Htmpx = zeros(length(f))
    RHtmpx = zeros(length(f))
    HRHtmpx = zeros(s.sv.n)

    b = zeros(s.sv.n)

    # x + W^½ * SB^½ * H' * (R \ (H * SB^½ * W^½ * x ))
    function fun!(x,fx)
        # tmpx = SB^½ W^½ x
        decompB!(s.sv,β,ivol,nus,nmax,α,sW * x,work1,work2,tmpx)

        # Htmpx = H * SB^½ W^½ x
        mul!(Htmpx,H,tmpx)

        @static if VERSION >= v"0.7.0-beta.0"
             # Htmpx = R \ (H * SB^½ W^½ x)
             ldiv!(R,Htmpx)

            # HRHtmpx = H' * (R \ (H * SB^½ W^½ x))
            mul!(HRHtmpx,H',Htmpx)
        else
            # Htmpx = R \ (H * SB^½ W^½ x)
            A_ldiv_B!(R,Htmpx)

            # HRHtmpx = H' * (R \ (H * SB^½ W^½ x))
            At_mul_B!(HRHtmpx,H,Htmpx)
        end

        # tmpx = SB^½ * H' * (R \ (H * SB^½ W^½ x ))
        decompB!(s.sv,β,ivol,nus,nmax,α,HRHtmpx,work1,work2,tmpx)

        # fx = x + W^½ SB^½ * H' * (R \ (H * SB^½ x ))
        mul!(fx,sW,tmpx)
        for i in 1:length(fx)
            fx[i] += x[i]
        end
    end

    # b = W^½ SB^{½} * H' * (R \ yo)
    decompB!(s.sv,β,ivol,nus,nmax,α,(H' * (R \ yo)),work1,work2,b)
    b = sW * b

    #@code_warntype decompB(n,s.sv,β,ivol,nus,nmax,α,H' * (R \ yo),work1,work2)

    # adjust tolerance
    tol = tol * s.sv.n / length(yo)

    # xp = (I + W^½ SB^½ * H' * (R^{-1} * (H * SB^½ * W^½)) )^{-1} b

    #@show DIVAnd.checksym(s.sv.n,fun!)

    xp,success,s.niter = DIVAnd.conjugategradient(
        fun!,b; tol = tol, maxit = maxit,
        progress = progress);

    # tmpx = SB^½ * W^½ * xp
    decompB!(s.sv,β,ivol,nus,nmax,α,sW * xp,work1,work2,tmpx)

    #@show mean(β),nmax,α

    # debugging function
    function decompB(x)

        #@show mean(β),nmax,α

        decompBx = similar(x)
        DIVAnd.decompB!(s.sv,β,ivol,nus,nmax,α,x,work1,work2,decompBx)

#        decompBx2 = (sW)^2 *  decompBx
#        DIVAnd.decompB!(s.sv,β,ivol,nus,nmax,α,decompBx2,work1,work2,decompBx)

        return decompBx
    end

    info = Dict(:B => x -> decompB((sW)^2 * decompB(x)),:H => H, :R => R, :yo => yo,
                :mask => mask, :ivol => ivol, :nus => nus, :β => β,
                :α => α, :nmax => nmax)
    return unpack(s.sv, tmpx,NaN)[1],s,info
end
