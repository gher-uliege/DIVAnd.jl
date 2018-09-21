if VERSION >= v"0.7.0-beta.0"
    using Test
else
    using Base.Test
end


#imax = 100
#jmax = 100

c = randn(10,10)
valex = -9999
c[3:5,6:10] = valex
cf = ufill(c,valex);
@test sum(cf == valex) == 0

c = randn(10,10,20)
valex = -9999
c[3:5,6:10,1:4] = valex
cf = ufill(c,valex);
@test sum(cf == valex) == 0


A = [1.,NaN,2.]
DIVAnd_fill!(A,NaN)
@test A ≈ [1.,1.5,2]


A = [1.,NaN,NaN,NaN,NaN]
DIVAnd_fill!(A,NaN)
@test A ≈ ones(size(A))
