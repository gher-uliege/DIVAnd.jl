# Return the analytical kernel and normalization factor.
#
# [mu,K] = divand_kernel(n,alpha)
# [mu,K] = divand_kernel(n,alpha,r)
#
# Analytical (normalized) kernels for infinite domain in dimension n and for
# coefficients alpha
# Input
#   n: number of dimensions
#   alpha: coefficients
#   r (optional): distance from origin
# Output:
#   K: kernel function evaluate at the values of r if present or a function handle
#   mu: normalization factor

function divand_kernel(n,alpha #,r
                       )
    # remove trailling zeros
    ind = maximum(find(!(alpha .== 0)))
    alpha = alpha[1:ind];

    m = length(alpha)-1;

    K = [];
    if [binomial(m,k) for k = 0:m] == alpha
        # alpha are binomial coefficients

        mu,K = divand_kernel_binom(n,m);
    else
        if [binomial(m,k) for k = 1:m] == alpha[2:ind]
            # alpha are binomial coefficients except first one
            warn("Semi-norm used?, check scaling $alpha")
            mu,K = divand_kernel_binom(n,m);
            #   correction for missing term CHECK IF NOT THE INVERSE
            #  Added fudge factor 2 to mimic same behaviour in test case
            jmscale=(1.0/2^(m))*sum(alpha[:])/2


            mu = mu*jmscale;

        else
            # unsupported sequence of alpha

            mu,K = divand_kernel_binom(n,m);
            warn("Unsupported norm used, check scaling $alpha")
            #   Scaling is correct of all alphas are binomials times a common factor

            jmscale=(1.0/2^(m))*sum(alpha[:])

            mu = mu*jmscale;
        end
    end

    # if nargin == 3
    #   # evaluate the kernel for the given values of r
    #   K = K(r);
    # end

    # if nargout == 3
    #   # determine at which distance r K(r) = 1/2
    #   rh = abs(fzero(@(r) K(abs(r))-.5,1));
    # end


    return mu,K #,rh

end



function divand_kernel_binom(n,m)


    # #if isequal(alpha,1)

    nu = m-n/2;
    mu = (4*pi)^(n/2) * gamma(m) / gamma(nu);
    K(x) = divand_rbesselk(nu,x);

    if nu <= 0
        warn("divand:nonorm","No normalization possible. Extend parameter alpha.");
        mu = 1;
    end

    return mu,K
end


function divand_rbesselk(nu,r)
    r = abs(r);

    if r == 0
        K = 1.
    else
        K = 2/gamma(nu) * ((r/2).^nu .* SpecialFunctions.besselk(nu,r));
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
