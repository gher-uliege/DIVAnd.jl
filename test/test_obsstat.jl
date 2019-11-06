using Test
import DIVAnd

nobs = 100
x = (randn(nobs), randn(nobs), randn(nobs))
v = randn(nobs)
ids = [randstring() for i = 1:nobs]

v[1] = NaN
v[2] = Inf

buf = IOBuffer()
@test_logs (:info, r".*Checking.*") match_mode = :any DIVAnd.checkobs(buf, x, v, ids)
output = lowercase(String(take!(buf)))

@test occursin("nan", output)
@test occursin("1", output)

@test occursin("inf", output)
@test occursin("2", output)


