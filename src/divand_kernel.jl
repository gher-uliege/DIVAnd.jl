"""
Return the analytical kernel and normalization factor.

mu,K = divand_kernel(n,alpha)

Analytical (normalized) kernels `K` for infinite domain in dimension `n` and for
coefficients `alpha` and normalization factor `mu`.
Input
  n: number of dimensions
  alpha: coefficients
Output:
  K(r): kernel function (function of the normalized distance `r`)
  mu: normalization factor
  len_scale: distance at which K(len_scale) = 0.6019072301972346 (which is besselk(1,1))
"""

function divand_kernel(n,alpha)
    # remove trailling zeros
    ind = maximum(find(!(alpha .== 0)))
    alpha = alpha[1:ind];

    m = length(alpha)-1;

    K = [];
    alpha_binomial = [binomial(m,k) for k = 0:m]

    if alpha_binomial == alpha
        # alpha are binomial coefficients

        mu,K = divand_kernel_binom(n,m);
    else
        if alpha_binomial[2:end] == alpha[2:end]
            # alpha are binomial coefficients except first one
            #warn("Semi-norm used?, check scaling $alpha")

            mu,K = divand_kernel_binom(n,m);
            #   correction for missing term CHECK IF NOT THE INVERSE
            #  Added fudge factor 2 to mimic same behaviour in test case
            jmscale=(1.0/2^(m))*sum(alpha[:])/2


            mu = mu*jmscale;

        else
            # unsupported sequence of alpha

            mu,K = divand_kernel_binom(n,m);
            #warn("Unsupported norm used, check scaling $alpha")
            #   Scaling is correct if all alphas are binomials times a common factor

            jmscale=(1.0/2^(m))*sum(alpha[:])

            mu = mu*jmscale;
        end
    end

    len_scale = Roots.fzero(x -> K(x) - SpecialFunctions.besselk(1,1),1)

    return mu,K,len_scale
end



function divand_kernel_binom(n,m)


    # #if isequal(alpha,1)

    nu = m-n/2;
    mu = (4*pi)^(n/2) * gamma(m) / gamma(nu);
    K(x) = divand_rbesselk(nu,x);

    if nu <= 0
        warn("divand:nonorm ","No normalization possible. Extend parameter alpha.");
        mu = 1;
    end

    return mu,K
end


function divand_rbesselk(nu,r)
    r = abs(r);

    if r == 0
        K = 1.
    else
        K = 2/gamma(nu) * ((r/2).^nu .* SpecialFunctions.besselk.(nu,r));
    end

    return K
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
