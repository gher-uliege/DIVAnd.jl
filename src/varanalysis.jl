
"""
Variational analysis similar to 3D-var

Kernel is the solution of the n-dimensional diffusion equation

∂c/∂t =  ∇ ⋅ (D ∇ c)

n-dimensional Green’s function

G(x,x',t) = (4πDt)^(-n/2)  exp( - |x -x'|² / (4Dt))
http://www.rpgroup.caltech.edu/~natsirt/aph162/diffusion_old.pdf

"""

function varanalysis(mask,pmn,xi,x,f,len,epsilon2; tol = 1e-5)
    n = ndims(mask)

    if isa(len,Number)
        len = ((fill(len,size(mask)) for i=1:n)...)
    elseif isa(len,Tuple)
        if isa(len[1],Number)
            len = ([fill(len[i],size(mask)) for i = 1:n]...)
        end

        for i=1:n
            if size(mask) != size(len[i])
                error("mask ($(formatsize(size(mask)))) and correlation length ($(formatsize(size(len[i])))) have incompatible size")
            end
        end
    end

    R = divand_obscovar(epsilon2,length(f));

    s = divand.divand_struct(mask)

    # observation constrain
    constrain = divand.divand_obs(s,xi,x,f,R)
    yo = constrain.yo
    H = constrain.H

    Ld = [mean(L) for L in len]
    nu = ([L.^2 for L in len]...)

    # D represents the Laplacian ∇ ⋅ (ν ∇ ϕ) where ν is the
    # correlation length-scale squared

    D = divand_laplacian(Val{:MatFun},mask,pmn,nu,falses(4))

    dx_min = 1/max([maximum(pm) for pm in pmn]...)
    nu_max = max([maximum(ν) for ν in nu]...)
    # need to proof this (currently this is just an analogy based on 2D)
    α = dx_min^2/(2 * nu_max * n)

    # 10% safety margin
    α = α / 1.1

    nmax = round(Int,1/(2*α))

    # the background error covariance matrix is
    # B =  (4π α nmax)^(n/2) prod(Ld) (I + α * D)^nmax;

    # x^T B x is the integral which takes also the volumn of each grid cell into account 

    sqrtivol = pack(s.sv, (sqrt(.*(pmn...)),))

    ivol,nus = divand_laplacian_prepare(mask,pmn,nu)

    function funB½(x)
        xup = unpack(s.sv,x)[1]

        for niter = 1:(nmax ÷ 2)
            xup += α * divand_laplacian_apply(ivol,nus,xup)
        end

        xup = (4π * α * nmax)^(n/4) * sqrt(prod(Ld)) * (sqrt.(ivol) .* xup)

        return pack(s.sv,(xup,))[:,1]
    end


    # matrix-like type
    B½ = MatFun(size(D),funB½,funB½)

    fun(x) = x + (B½ * (H' * (R \ (H * (B½ * x)))))
    b = B½ * (H' * (R \ yo))

    # adjust tolerance
    tol = tol * s.sv.n / length(yo)

    xp,success,niter = divand.conjugategradient(fun,b; tol = tol);
    xa = B½ * xp

    return unpack(s.sv,xa)[1]
end
