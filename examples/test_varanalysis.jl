# A simple example of divand in 4 dimensions
# with observations from an analytical function.

using divand
using PyPlot

# observations
xy = ([.5],[.5])
f = [1]

# final grid
testsizexy=11
testsizez=11
testsizet=11

testsizexy=21
testsizexy=31
testsizexy=41
#testsizexy=101
testsizez=11
testsizet=11

# mask: all points are valid points
# this problem has a simple cartesian metric
# pm is the inverse of the resolution along the 1st dimension,...
mask,pmn,xyi = divand_rectdom(linspace(0,1,testsizexy),
                              linspace(0,1,testsizexy))



sv = statevector((mask,))

# correlation length
lenxy = (0.2,0.2)

# obs. error variance normalized by the background error variance
epsilon2 = 1;

# fi is the interpolated field

#@time fi,sd = divandrun(mask,pmn,xyi,xy,f,(lenx,leny,0,0),epsilon2);

#@time fip,sp = divandrun(mask3,pmn,xyi,(lon,lat,depth,time2),va,(1,1,0,0.00001),epsilon2)
#fip = fip+vm;

# tolerance on the gradient A x - b
tol = 1e-1
tol = 0.5
# tolerance on the result x
tolres = 1e-3

kwargs = [(:tol, tol),(:maxit,10000),(:minit,0)]
kwargs = [(:tol, tol),(:maxit,100),(:minit,0)]
kwargs = [(:tol, tol),(:maxit,10),(:minit,0)]

function progressiter(iter,x,r,tol2,fun,b)
    # J = x' A x/2 - b'*x  
    
    J = x ⋅ fun(x)/2 - b ⋅ x

    @show iter, r ⋅ r, tol2, minimum(x),maximum(x),J

    return nothing
end

function progress(iter,x,r,tol2,fun,b)
    
    progressiter(iter,x,r,tol2,fun,b)
#    return nothing

    f, = unpack(sv,x)
    fr, = unpack(sv,r)
    #m = (size(fi2,1) + 1) ÷ 2
    m = ([size(mask)...] + 1) ÷ 2

    ca = [-0.1,.5]
    clf()
    subplot(1,2,1); 
    title("iteration $(iter)")
    pcolor(f); xlabel("x"), ylabel("y");
    clim(ca)
    colorbar()

    subplot(1,2,2); 
    title("iteration $(iter) residual")
    pcolor(fr); xlabel("x"), ylabel("y");
    colorbar()
    
    savefig(@sprintf("test2D-frame%04d.png",iter))

    iter = iter + 1
    
    @show iter
end

compPC2D(iB,H,R) = x -> s.P*x

#D = divand_laplacian(Val{:sparse},mask,pmn,lenxy,falses(4))


nu = ([L.^2 for L in lenxy]...)

#D = divand_laplacian(Val{:MatFun},mask,pmn,nu,falses(4))
D = divand_laplacian(Val{:sparse},mask,pmn,nu,falses(4))


# http://www.rpgroup.caltech.edu/~natsirt/aph162/diffusion_old.pdf
# equation 60
# n-dimensional Green’s function
# ∂c  
# --  =  ∇ ⋅ (D ∇ c)
# ∂t 

# G(x,x',t) = (4πDt)^(-n/2)  exp( - |x -x'|² / (4Dt))
# D = L²

#    The function `compPC` returns the
#    preconditioner `fun(x)` representing `M \ x` (the inverse of M times x)
#    where `M` is a positive defined symmetric matrix [1].

#    Effectively, the system E⁻¹ A (E⁻¹)ᵀ (E x) = E⁻¹ b is solved for (E x) where E Eᵀ = M.
#    Ideally, M should this be similar to A, so that E⁻¹ A (E⁻¹)ᵀ is close to the identity matrix.

# Here
#  A x = b
#  P
# need to match diva kernel to Gaussian

function compPC(iB,H,R) 

    # 4 L² dt * nmax ~ 2 L²
    # nmax ~  1/(2 dt)
    dt = 0.005
    dt = 0.001
    # need to tune the factor
    nmax = round(Int,1/dt/2)
   
    @show nmax
    nmax = 26
    nmax = 10
    nmax = 30
    nmax = 5
    function pc(x)        
        for i = 1:nmax
            x += dt * D*x
        end
        return x
    end

    return pc
end


@time fid,sd = divandrun(mask,pmn,xyi,xy,f,lenxy,epsilon2)

# @time fi2,s = divandrun(mask,pmn,xyi,xy,f,lenxy,epsilon2;
#                         kwargs...,inversion=:pcg,operatortype=Val{:sparse}
# #                        ,fi0=fi
#                         ,compPC = compPC
#                         ,progress = progress
# #                        ,progress = progressiter
#                          )
# @show s.niter

fullP = sd.P * eye(size(sd.P,1))

nmax = 100;
α = 0.005
B = (eye(size(D,1)) + α * D)^nmax;
B½ = (eye(size(D,1)) + α * D)^(nmax ÷ 2)

H = sd.obsconstrain.H
R = sd.obsconstrain.R
yo = sd.obsconstrain.yo

iR = Diagonal(1./diag(R))
invPa = Hermitian(inv(B) + H' * (R \ H))
xa = invPa \ (H' * (R \ yo))


I = eye(size(B,1))
Pp = I + B½ * H' * (R \ H) * B½
xp = Pp \ (B½ * (H' * (R \ yo)))
xa2 = B½ * xp;

fun(x) = x + (B½ * (H' * (R \ (H * (B½ * x)))))
b = B½ * (H' * (R \ yo))
xp,success,niter = divand.conjugategradient(fun,b,tol = 1e-5);
xa3 = B½ * xp




function varanalysis(mask,pmn,xi,x,f,lenxy,epsilon2; tol = 1e-5)

    R = divand_obscovar(epsilon2,length(f));

    s = divand.divand_struct(mask)
    # observation constrain

    @show size(xi[1])
    @show size(s.mask)


    constrain = divand.divand_obs(s,xi,x,f,R)
    yo = constrain.yo
    H = constrain.H

    nu = ([L.^2 for L in lenxy]...)

    # D represents the Laplacian ∇ ⋅ (ν ∇ ϕ) where ν is the 
    # correlation length-scale squared

    D = divand_laplacian(Val{:MatFun},mask,pmn,nu,falses(4))

    α = 0.001
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
    xp,success,niter = divand.conjugategradient(fun,b; tol = tol);
    xa = B½ * xp

    return xa
end

@time xa4 = varanalysis(mask,pmn,xyi,xy,f,lenxy,epsilon2; tol = 1e-5)


# Copyright (C) 2014, 2017 Alexander Barth <a.barth@ulg.ac.be>
#                          Jean-Marie Beckers <JM.Beckers@ulg.ac.be>
#
# This program is free software; you can redistribute it and/or modify it under
# the terms of the GNU General Public License as published by the Free Software
# Foundation; either version 2 of the License, or (at your option) any later
# version.
#
# This program is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
# FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
# details.
#
# You should have received a copy of the GNU General Public License along with
# this program; if not, see <http://www.gnu.org/licenses/>.
