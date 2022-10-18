using Test
using Random
import DIVAnd
using StableRNGs

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


lon = [1,2,1]
lat = [10,20,10]
val = [1,2,-1]
ulon,ulat = DIVAnd.statpos((lon, lat))
@test sort(ulon) ≈ [1,2]
@test sort(ulat) ≈ [10,20]

(ulon,ulat),meanval,stdval,count = DIVAnd.statpos(val, (lon, lat))
@test sort(meanval) ≈ [0, 2]

# test DIVAnd.randsplit

rng = StableRNG(123)
nobs = 100
x = (rand(rng,1:5,nobs), rand(rng,1:5,nobs))

fractions = (0.9, 0.1)

groupindex = DIVAnd.randsplit(x,fractions)

# check that no observations with the same coordinates is split accross
# different groups.

groupindex_map = zeros(Int,5,5)
for j = 1:5
    for i = 1:5
        for l = 1:length(x[1])
            if (i == x[1][l]) && (j == x[2][l])
                if groupindex_map[i,j] == 0
                    groupindex_map[i,j] = groupindex[l]
                else
                    # location already associated to a group,
                    # check that this association is consistent
                    #@show groupindex_map[i,j], groupindex[l]
                    @assert groupindex_map[i,j] == groupindex[l]
                end
            end
        end
    end
end

