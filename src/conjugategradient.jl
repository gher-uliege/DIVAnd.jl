
"""
    pc_none!(x,fx)

Dummy call-back function when no preconditioner is used. `fx` will be equal to `x`.
"""
function pc_none!(x, fx)
    fx[:] = x
end

"""
    xAy, yATx = checksym(n,fun!)

Check if the the function `fun!` represents a symmetric matrix when applied on
random vectors of size `n`.
"""
function checksym(n, fun!)
    x = randn(n)
    y = randn(n)
    Ax = zeros(n)
    Ay = zeros(n)

    fun!(x, Ax)
    fun!(y, Ay)

    return (y ⋅ Ax), (x ⋅ Ay)
end

function cgprogress(iter, x, r, tol2, fun!, b)
    if iter == 1
        print("|  Iteration | Cost function | Norm of residual/√n | max. norm/√n |\n")
        print("|------------|---------------|---------------------|--------------|\n")
    end

    # this is the same
    # Ax = zeros(size(x))
    # fun!(x,Ax)

    Ax = b - r
    J = (x ⋅ Ax) / 2 - (b ⋅ x)
    n = length(r)

    print("|")
    printstyled("$(@sprintf("%11d",iter))", bold = true)
    print(" |")
    printstyled("$(@sprintf("%14.3f",J))", color = :light_magenta)
    print(" |")
    printstyled("$(@sprintf("%20f",sqrt((r ⋅ r)/n)))", color = :red)
    print(" |")
    printstyled("$(@sprintf("%13f",sqrt(tol2/n)))", color = :default)
    print(" |\n")
end


"""
    x,cgsuccess,niter = conjugategradient(fun!,b)

Solve a linear system with the preconditioned conjugated-gradient method:
A x = b
where `A` is a symmetric positive defined matrix and `b` is a vector.
Equivalently the solution `x` minimizes the cost function
J(x) = ½ xᵀ A x - bᵀ x.

The function `fun!(x,fx)` computes fx which is equal to  `A*x`.
For example:

```
function fun!(x,fx)
    fx[:] = A*x
end
```

Note that the following code will NOT work, because a new array `fx` would be created and it would not be passed back to the caller.

```
function fun!(x,fx)
    fx = A*x # bug!
end
```
The function `fun!` works in-place to reduce the amount of memory allocations.

# Optional input arguments
* `x0`: starting vector for the interations
* `tol`: tolerance on  |Ax-b| / |b|
* `maxit`: maximum of interations
* `pc!`: the preconditioner. The functions `pc(x,fx)` computes fx = M⁻¹ x (the inverse of M times x) where `M` is a symmetric positive defined matrix. Effectively, the system E⁻¹ A (E⁻¹)ᵀ (E x) = E⁻¹ b is solved for (E x) where E Eᵀ = M. Ideally, M should this be similar to A, so that E⁻¹ A (E⁻¹)ᵀ is close to the identity matrix. The function `pc!` should be implemented in a similar way than `fun!` (see above).

# Output
* `x`: the solution
* `cgsuccess`: true if the interation converged (otherwise false)
* `niter`: the number of iterations
"""
function conjugategradient(
    fun!,
    b::Vector{T};
    x0::Vector{T} = zeros(size(b)),
    tol::T = 1e-6,
    maxit::Int = min(size(b, 1), 100000),
    minit::Int = 0,
    pc! = pc_none!,
    progress = (iter, x, r, tol2, fun!, b) -> nothing,
) where {T}

    success = false
    n = length(b)

    bb = b ⋅ b

    if bb == 0
        GC.enable(true)
        return zeros(size(b)), true, 0
    end

    # relative tolerance
    tol2 = tol^2

    # absolute tolerance
    tol2 = tol2 * bb

    # initial guess
    x = x0

    # memory allocation
    Ap = similar(x)
    z = similar(x)

    # gradient at initial guess
    fun!(x, Ap)
    r = b - Ap

    # quick exit
    if r ⋅ r < tol2
        GC.enable(true)
        return x, true, 0
    end

    # apply preconditioner
    pc!(r, z)

    # first search direction == gradient
    p = copy(z)

    # compute: r' * inv(M) * z (we will need this product at several
    # occasions)

    # ⋅ is the dot vector product and returns a scalar
    zr_old = r ⋅ z

    alpha = zeros(T, maxit)
    beta = zeros(T, maxit + 1)

    kfinal = maxit
    #gc()
    for k = 1:maxit
        # compute A*p
          #@show k
        fun!(p, Ap)

        # how far do we need to go in direction p?
        # alpha is determined by linesearch
        alpha[k] = zr_old / (p ⋅ Ap)

        # get new estimate of x
        # x = x + alpha[k]*p
        x = BLAS.axpy!(alpha[k], p, x)

        # recompute gradient at new x. Could be done by
        # r = b-fun(x)
        # but this does require an new call to fun
        # r = r - alpha[k]*Ap
        r = BLAS.axpy!(-alpha[k], Ap, r)

        progress(k, x, r, tol2, fun!, b)

        #if mod(k,20)==1
        #    @show k, r ⋅ r,tol2,size(r)
        #end

        if ((r ⋅ r) < tol2) && (k >= minit)
            success = true
            #@show k
            kfinal = k
            break
        end

        # apply pre-conditionner

        pc!(r, z)

        zr_new = r ⋅ z

        # Fletcher-Reeves
        beta[k+1] = zr_new / zr_old
        # Polak-Ribiere
        # beta[k+1] = r'*(r-r_old) / zr_old
        # Hestenes-Stiefel
        # beta[k+1] = r'*(r-r_old) / (p'*(r-r_old))
        # beta[k+1] = r'*(r-r_old) / (r_old'*r_old)

        # p = z + beta[k+1]*p
        for i = 1:n
            p[i] = z[i] + beta[k+1] * p[i]
        end

        zr_old = zr_new
    end

    GC.enable(true)

    return x, success, kfinal

end

# Copyright (C) 2004,2017  Alexander Barth           <a.barth@ulg.ac.be>
#                          Jean-Marie Beckers          <jm.beckers@ulg.ac.be>
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
