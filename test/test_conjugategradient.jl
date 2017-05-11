using Base.Test 

n = 10; # dimension

S = randn(n,n);
b = randn(n);

A = 10*eye(n,n) + 0.01*S*S'; # symmetric and positive defined matrix
function fun!(x,fx)
    fx[:] = A*x
end

function pc_exact!(x,fx)
    fx[:] = A\x
end

function pc_jacobi!(x,fx)
    fx[:] = Diagonal(diag(A))\x
end


tol = 1e-4
kwargs = [(:tol, tol)]

x,success,niter0 = divand.conjugategradient(fun!,zeros(b); kwargs...)
@test x ≈ zeros(b)

x,success,niter1 = divand.conjugategradient(fun!,b; kwargs...)
@test norm(A*x - b)/norm(b) < tol

x,success,niter2 = divand.conjugategradient(fun!,b; kwargs..., pc! = pc_exact!)
@test norm(A*x - b)/norm(b) < tol
@test niter2 == 1

x,success,niter3 = divand.conjugategradient(fun!,b; kwargs..., pc! = pc_jacobi!)
@test norm(A*x - b)/norm(b) < tol
@test niter3 <= niter1



