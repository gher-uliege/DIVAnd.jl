# Testing DIVAnd in 2 dimensions

if VERSION >= v"0.7.0-beta.0"
    using Test
else
    using Base.Test
end

# grid of background field
mask,(pm,pn),(xi,yi) = DIVAnd_squaredom(
    2,range(0, stop = 1, length = 15))

epsilon = 1e-10;

# grid of observations
x,y = ndgrid(
    range(epsilon, stop = 1-epsilon, length = 5),
    range(epsilon, stop = 1-epsilon, length = 5))

x = x[:]
y = y[:]
v = sin.(x*6) .* cos.(y*6)

# correlation length
lenx = .15;
leny = .15;

epsilon2 = 1.

# tolerance on the gradient A x - b
tol = 1e-4
# tolerance on the result x
tolres = 1e-3

kwargs = [(:tol, tol),(:maxit,1000)]

# Try different solvers
# * compare results to direct solver
# * verify that preconditioners reduce the number of iterations

# direct inversion
va_chol,s = DIVAndrun(mask,(pm,pn),(xi,yi),(x,y),v,(lenx,leny),epsilon2;
                      kwargs..., inversion=:chol)

for jj=1:4
    # iterative (without preconditioner)
    va_iter,s_np = DIVAndrun(mask,(pm,pn),(xi,yi),(x,y),v,(lenx,leny),epsilon2;
                             kwargs..., inversion=:pcg,btrunc=jj)
    @test norm(va_chol[mask] - va_iter[mask]) < tolres
end

# iterative (without preconditioner)
va_iter,s_np = DIVAndrun(mask,(pm,pn),(xi,yi),(x,y),v,(lenx,leny),epsilon2;
                         kwargs..., inversion=:pcg)
@test norm(va_chol[mask] - va_iter[mask]) < tolres

# iterative (without preconditioner but a good starting point :-)
va_iter,s = DIVAndrun(mask,(pm,pn),(xi,yi),(x,y),v,(lenx,leny),epsilon2;
                      kwargs..., inversion=:pcg, fi0 = va_chol)
@test va_iter ≈ va_chol
@test s.niter == 0


# iterative (with preconditioner)

va_iter,s_wp = DIVAndrun(mask,(pm,pn),(xi,yi),(x,y),v,(lenx,leny),epsilon2;
                         kwargs..., inversion=:pcg,compPC = DIVAnd_pc_sqrtiB)
@test norm(va_chol[mask] - va_iter[mask]) < tolres
@test s_wp.niter < s_np.niter


# iterative (with custom preconditioner)
function compPC(iB,H,R)
    function fun!(x,fx)
        fx[:] = iB \ x
    end
    return fun!
end
va_iter,s_wp = DIVAndrun(mask,(pm,pn),(xi,yi),(x,y),v,(lenx,leny),epsilon2;
                         kwargs..., inversion=:pcg,compPC = compPC);
@test norm(va_chol[mask] - va_iter[mask]) < tolres

# iterative dual (without precondiditioner)

va_dual,s_np = DIVAndrun(mask,(pm,pn),(xi,yi),(x,y),v,(lenx,leny),epsilon2;
                         kwargs..., primal=false);
@test norm(va_chol[mask] - va_dual[mask]) < tolres

# iterative dual (with precondiditioner)
# This is not efficient for large cases, only a consistency check
function compPCdual(iB,H,R)
    B = CovarIS(iB)
    factorize!(B)

    M = H * (B * H') + sparse_diag(diag(R));

    iM = CovarIS(M)
    factorize!(iM)

    function fun!(x,fx)
        fx[:] = iM * x
    end
    return fun!
end
va_dual,s_wp = DIVAndrun(mask,(pm,pn),(xi,yi),(x,y),v,(lenx,leny),epsilon2;
                         kwargs..., primal=false,compPC = compPCdual);
@test norm(va_chol[mask] - va_dual[mask]) < tolres
@test s_wp.niter < s_np.niter



# Copyright (C) 2014, 2016 Alexander Barth <a.barth@ulg.ac.be>
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
