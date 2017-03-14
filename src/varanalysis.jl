"""
Variational analysis similar to 3D-var

Kernel is the solution of the n-dimensional diffusion equation

∂c/∂t = =  ∇ ⋅ (D ∇ c)
 
n-dimensional Green’s function

G(x,x',t) = (4πDt)^(-n/2)  exp( - |x -x'|² / (4Dt))

"""

function varanalysis(mask,pmn,xi,x,f,lenxy,epsilon2; tol = 1e-5)

    R = divand_obscovar(epsilon2,length(f));

    s = divand.divand_struct(mask)

    # observation constrain
    constrain = divand.divand_obs(s,xi,x,f,R)
    yo = constrain.yo
    H = constrain.H

    nu = ([L.^2 for L in lenxy]...)

    # D represents the Laplacian ∇ ⋅ (ν ∇ ϕ) where ν is the 
    # correlation length-scale squared

    D = divand_laplacian(Val{:MatFun},mask,pmn,nu,falses(4))

    dx_min = 1/max([maximum(pm) for pm in pmn]...)
    nu_max = max([maximum(ν) for ν in nu]...)
    # need to proof this (currently this is just an analogy based on 2D)
    α = dx_min^2/(2 * nu_max * s.n)

    # 10% safety margin
    α = α / 1.1

    nmax = round(Int,1/(2*α))

    # the background error covariance matrix is    
    # B = (I + α * D)^nmax;

    # the square root of the background error covariance matrix:
    # B^{1/2} = (I + α * D)^(nmax/2);


    function funB½(x)
        for n = 1:(nmax ÷ 2)
            x += α * (D*x)
        end
        return x
    end

    # matrix-like type of 
    B½ = MatFun(size(D),funB½,funB½)

    fun(x) = x + (B½ * (H' * (R \ (H * (B½ * x)))))
    b = B½ * (H' * (R \ yo))
    #@show maximum(b)

    xp,success,niter = divand.conjugategradient(fun,b; tol = tol);
    xa = B½ * xp

    return unpack(s.sv,xa)[1]
end
