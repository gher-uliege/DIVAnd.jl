using Base.Test 

n = 10; # dimension

S = randn(n,n);
b = randn(n);

A = 10*eye(n,n) + 0.01*S*S'; # symmetric and positive defined matrix
fun(x) = A*x

tol = 1e-4
kwargs = [(:tol, tol)]

x,success,niter0 = divand.conjugategradient(fun,zeros(b); kwargs...)
@test x ≈ zeros(b)

x,success,niter1 = divand.conjugategradient(fun,b; kwargs...)
@test norm(A*x - b)/norm(b) < tol

x,success,niter2 = divand.conjugategradient(fun,b; kwargs..., pc = x -> A\x)
@test norm(A*x - b)/norm(b) < tol
@test niter2 == 1

x,success,niter3 = divand.conjugategradient(fun,b; kwargs..., pc = x -> Diagonal(diag(A))\x)
@test norm(A*x - b)/norm(b) < tol
@test niter3 <= niter1



