"""
    mu,K,len_scale = DIVAnd_kernel(n,alpha)

Return the analytical kernel and normalization factor.

Analytical (normalized) kernels `K` for infinite domain in dimension `n` and for
coefficients `alpha` and normalization factor `mu`.

`K(r)` is the kernel function (function of the normalized distance `r`),
`len_scale` is the distance at which `K(len_scale)` = 0.6019072301972346 (which is besselk(1,1))
"""
function DIVAnd_kernel(n,alpha)
    # precision on len_scale
    eps = 1e-8

    # remove trailling zeros
    ind = findlast(alpha .!= 0)
    alpha = alpha[1:ind]

    m = length(alpha)-1;

    alpha_binomial = [binomial(m,k) for k = 0:m]

    mu,K = DIVAnd_kernel_binom(n,m);

    if alpha_binomial == alpha
        # alpha are binomial coefficients
        # ok pass
    else
        if alpha_binomial[2:end] == alpha[2:end]
            # alpha are binomial coefficients except first one
            #@warn "Semi-norm used?, check scaling $alpha"

            #   correction for missing term CHECK IF NOT THE INVERSE
            #  Added fudge factor 2 to mimic same behaviour in test case
            jmscale=(1.0/2^(m))*sum(alpha[:])/2


            mu = mu*jmscale;

        else
            # unsupported sequence of alpha

            #@warn "Unsupported norm used, check scaling $alpha"
            #   Scaling is correct if all alphas are binomials times a common factor

            jmscale=(1.0/2^(m))*sum(alpha[:])

            mu = mu*jmscale;
        end
    end

    len_scale,maxiter = fzero(x -> K(x) - SpecialFunctions.besselk(1,1),0.,100.,eps)

    return mu,K,len_scale
end



function DIVAnd_kernel_binom(n,m)
    nu = m-n/2
    mu = (4*pi)^(n/2) * gamma(m) / gamma(nu)
    K(x) = DIVAnd_rbesselk(nu,x)

    if nu <= 0
        @warn "DIVAnd:nonorm ","No normalization possible. Extend parameter alpha."
        mu = 1.
    end

    return mu,K
end


function DIVAnd_rbesselk(nu,r)
    r = abs(r);

    if r == 0
        K = 1.
    else
        K = 2/gamma(nu) * ((r/2).^nu .* SpecialFunctions.besselk.(nu,r));
    end

    return K
end


"""
fzero(f,x0,x1,eps; maxiter = Inf)
find the zero of the function f between x0 and x1 assuming x0 < x1 at a precision eps.
"""
function fzero(f,x0,x1,eps; maxiter = 1000)

    if x0 > x1
        x0,x1 = x1,x0
    end

    fx0 = f(x0);
    fx1 = f(x1);
    xc = (x0+x1)/2.;
    fxc = f(xc);
    niter = 0

    if fx0*fx1 > 0
        error("function at x0 (f($x0)=$(fx0)) and x1 (f($x1)=$(fx1)) should have a different sign")
    end

    while (x1-x0 > eps) && fxc != 0 && niter < maxiter
        if fx0*fxc > 0
            x0 = xc
            fx0 = fxc
        else
            x1 = xc
            fx1 = fxc
        end

        xc = (x0+x1)/2.
        fxc = f(xc)
        niter += 1
    end

    return xc,niter
end

# Copyright (C) 2014, 2017 Alexander Barth <a.barth@ulg.ac.be>
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
