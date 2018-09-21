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
B = similar(A)
DIVAnd_fill!(A,B,NaN)
@test B ≈ [1.,1.5,2]
