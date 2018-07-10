
fun(x) = x^2-2
x0,maxiter = DIVAnd.fzero(fun, 0.,10.,1e-8)

@test abs(fun(x0)) < 1e-5
@test_throws ErrorException DIVAnd.fzero(fun, 0.,0.1,1e-5)
