if VERSION >= v"0.7.0-beta.0"
    using Test
else
    using Base.Test
end
import DIVAnd

nobs = 100
x = (randn(nobs),randn(nobs),randn(nobs))
v = randn(nobs)
ids = [randstring() for i in 1:nobs]

v[1] = NaN
v[2] = Inf

buf = IOBuffer()
@static if VERSION >= v"0.7.0"
    @test_logs (:info,r".*Checking.*") match_mode=:any DIVAnd.checkobs(buf,x,v,ids)
else
    @test_warn r".*Checking.*" DIVAnd.checkobs(buf,x,v,ids)
end
output = lowercase(String(take!(buf)))

@test occursin("nan",output)
@test occursin("1",output)

@test occursin("inf",output)
@test occursin("2",output)


