using Base.Test
import DIVAnd

nobs = 100
x = (randn(nobs),randn(nobs),randn(nobs))
v = randn(nobs)
ids = [randstring() for i in 1:nobs]
 
v[1] = NaN
v[2] = Inf

buf = IOBuffer()
DIVAnd.checkobs(buf,x,v,ids)

output = lowercase(String(take!(buf)))

@test contains(output,"nan")
@test contains(output,"1")

@test contains(output,"inf")
@test contains(output,"2")


