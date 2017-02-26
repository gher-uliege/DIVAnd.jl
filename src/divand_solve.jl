"""
Solve the variational problem.

fi = divand_solve(s)

Derive the analysis based on all contraints included in s and using the
observations yo

Input:
  s: structure created by divand_factorize
  fi0: starting point for iterative primal methods
  f0: starting point for the iterative dual method

Output:
  fi: analyzed field
"""

function divand_solve!(s,fi0,f0)

    H = s.H;
    sv = s.sv;
    R = s.R;
    yo = s.yo;


    if s.primal
        if s.inversion == :chol
            P = s.P;
            fpi =  P * (H'* (R \ yo[:]));
        else
            HiRyo = H'* (R \ yo[:]);
            fun(x) = s.iB*x + H'*(R \ (H * x));

            fpi,success,s.niter = conjugategradient(fun,HiRyo,tol = s.tol,
                                                    maxit = s.maxit,
                                                    minit = s.minit,
                                                    x0 = fi0,
                                                    pc = s.preconditioner,
                                                    progress = s.progress
                                                    )

            if !success
                warn("Preconditioned conjugate gradients method did not converge")
            end

            #s.P = CovarLanczos(Q,T);
        end
    else # dual
        B = CovarIS(s.iB);

        # fun(x) computes (H B H' + R)*x
        fundual(x) = H * (B * (H'*x)) + R*x;

        tmp,success,s.niter = conjugategradient(fundual,yo,tol = s.tol,
                                                maxit = s.maxit,
                                                minit = s.minit,
                                                x0 = f0,
                                                pc = s.preconditioner,
                                                progress = s.progress
                                                )
        if !success
            warn("Preconditioned conjugate gradients method did not converge")
        end

        fpi = B * (H'*tmp);
    end


    fi, = statevector_unpack(sv,fpi);

    fi[!s.mask] = NaN;

    return fi
end

# Copyright (C) 2014,2017 Alexander Barth <a.barth@ulg.ac.be>
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
