# Testing divand in 2 dimensions

using Base.Test

# grid of background field
xi,yi = ndgrid(linspace(0,1,30),linspace(0,1,30))

mask = trues(size(xi))
pm = ones(size(xi)) / (xi[2,1]-xi[1,1])
pn = ones(size(xi)) / (yi[1,2]-yi[1,1])

epsilon = 1e-10;

# grid of observations
x,y = ndgrid(linspace(epsilon,1-epsilon,10),linspace(epsilon,1-epsilon,10))
x = x[:]
y = y[:]
v = sin(x*6) .* cos(y*6)

# correlation length
lenx = .15;
leny = .15;

epsilon2 = 1;

# tolerance on the gradient A x - b
tol = 1e-4
# tolerance on the result x
tolres = 1e-3

kwargs = [(:tol, tol),(:maxit,1000)]

# Try different solvers
# * compare results to direct solver
# * verify that preconditioners reduce the number of iterations

# direct inversion
va_chol,s = divandrun(mask,(pm,pn),(xi,yi),(x,y),v,(lenx,leny),epsilon2;
                            kwargs..., inversion=:chol)

# iterative (without preconditioner)
va_iter,s_np = divandrun(mask,(pm,pn),(xi,yi),(x,y),v,(lenx,leny),epsilon2;
                            kwargs..., inversion=:pcg)
@test norm(va_chol[mask] - va_iter[mask]) < tolres

# iterative (without preconditioner but a good starting point :-)
va_iter,s = divandrun(mask,(pm,pn),(xi,yi),(x,y),v,(lenx,leny),epsilon2;
                            kwargs..., inversion=:pcg, fi0 = va_chol)
@test va_iter ≈ va_chol
@test s.niter == 0


# iterative (with preconditioner)

va_iter,s_wp = divandrun(mask,(pm,pn),(xi,yi),(x,y),v,(lenx,leny),epsilon2;
                      kwargs..., inversion=:pcg,compPC = divand_pc_sqrtiB)
@test norm(va_chol[mask] - va_iter[mask]) < tolres
@test s_wp.niter < s_np.niter


# iterative (with custom preconditioner)
function compPC(iB,H,R)
    return x -> iB \ x;
end
va_iter,s_wp = divandrun(mask,(pm,pn),(xi,yi),(x,y),v,(lenx,leny),epsilon2;
                      kwargs..., inversion=:pcg,compPC = compPC);
@test norm(va_chol[mask] - va_iter[mask]) < tolres

# iterative dual (without precondiditioner)

va_dual,s_np = divandrun(mask,(pm,pn),(xi,yi),(x,y),v,(lenx,leny),epsilon2;
                      kwargs..., primal=false);
@test norm(va_chol[mask] - va_dual[mask]) < tolres

# iterative dual (with precondiditioner)
function compPCdual(iB,H,R)
    B = CovarIS(iB)
    M = H * (B * H') + sparse_diag(diag(R));
    iM = CovarIS(M);
    factorize!(iM);
    return x -> iM * x;
end
va_dual,s_wp = divandrun(mask,(pm,pn),(xi,yi),(x,y),v,(lenx,leny),epsilon2;
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
