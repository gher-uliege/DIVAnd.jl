


"""
    ulon,ulat = statpos(lon,lat)
    ulon,ulat = statpos((lon,lat,...))

Return unique positions (`ulon`, `ulat`) as well their mean,
standard deviation and count of the vector of observations `val` located
at the positions `lon` and `lat`.

"""
function statpos(x::NTuple)
    allpos = collect(zip(x...))
    uniquepos = collect(Set(allpos))
    return ntuple(i -> getindex.(uniquepos,i),length(x))
end

@deprecate statpos(lon::AbstractVector, lat::AbstractVector) statpos((lon,lat))



"""
    (ulon,ulat),meanval,stdval,count = statpos(val,(lon,lat))
    (ulon,ulat,...),meanval,stdval,count = statpos(val,(lon,lat,...))

Return unique positions (`ulon`, `ulat`) as well as their mean,
standard deviation and count of the vector of observations `val` located
at the positions `lon` and `lat`.

"""
function statpos(val, x::NTuple)
    allpos = collect(zip(x...))
    uniquepos = collect(Set(allpos))

    count = zeros(Int, length(uniquepos))
    sumval = zeros(length(uniquepos))
    sumval2 = zeros(length(uniquepos))

    mapping = Dict([(up,i) for (i,up) in enumerate(uniquepos)])
    for j = 1:length(allpos)
        i = mapping[allpos[j]]
        count[i] += 1
        sumval[i] += val[j]
        sumval2[i] += val[j]^2
    end

    meanval = sumval ./ count
    stdval = sumval2 ./ count - meanval .^ 2
    stdval[stdval.<0] .= 0
    stdval = sqrt.(stdval)

    uniquex = ntuple(i -> getindex.(uniquepos,i),length(x))

    return uniquex, meanval, stdval, count
end

@deprecate statpos(val, lon::AbstractVector, lat::AbstractVector) begin
    (ulon,ulat),meanval,stdval,count = statpos(val,(lon,lat))
    return ulon, ulat, meanval, stdval, count
end

"""
    checkobs(x,v,ids)
    checkobs(io::IO,x,v,ids)

Print some basic information about the coordinates `x` (tuple of vector) and
values `v` (vector) having the identifier `ids` (vector of strings) to check
erroneous data. It prints wheter NaNs or Infs are found and the minimum and
maximum value.

If the argument `io` is provided, the information is input/output stream `io`.
"""
checkobs(x, v, ids) = checkobs(stdout, x, v, ids)

function checkobs(io::IO, x, v, ids)
    @info "Checking ranges for dimensions and observations"
    function check(xc, ids, name)
        fmt = "%55s"

        if eltype(xc) <: AbstractFloat

            for (fun, str) in [(isnan, "NaN"), (isinf, "infinity")]
                n = sum(fun, xc)

                if n > 0
                    @printf(io, "%55s", "number of values equal to $(str) in $(name): ")
                    printstyled(io, n, color = :red)

                    j = findfirst(isnan, xc)
                    printstyled(
                        io,
                        " [first value at index $(j) and id $(ids[j])]",
                        bold = true,
                    )
                    println(io)
                end
            end
        end

        @printf(io, "%55s", "minimum and maximum of $(name): ")
        print(io, extrema(xc[isfinite.(xc)]))
        println(io)

    end
    # loop over all dimensions
    for i = 1:length(x)
        check(x[i], ids, "obs. dimension $(i)")
    end

    check(v, ids, "data")
end



"""
    groupindex = DIVAnd.randsplit(x,fractions)

Split the observations based on their positions `x` (a tuple of vectors with
the coordinates) into as many groups as their are elements in `fractions`
(which are the approximate fraction of unique positions for the different group).
Observations with the same coordinates will not be split accross different groups.

### Example

In this example, split the data such that 90 % belong to the analysis dataset
and 10% to validation dataset.

```julia
fractions = (0.9,0.1)
# random data
lonobs = rand(1:5,100)
latobs = rand(1:5,100)
value = rand(100)
groupindex = DIVAnd.randsplit((lonobs,latobs),fractions)
lonobs_analysis = lonobs[groupindex .== 1]
latobs_analysis = latobs[groupindex .== 1]
value_analysis = value[groupindex .== 1]

lonobs_validation = lonobs[groupindex .== 2]
latobs_validation = latobs[groupindex .== 2]
value_validation = value[groupindex .== 2]
```
"""
function randsplit(x,fractions; rng = Random.GLOBAL_RNG)
    # cumulative fractions
    @static if VERSION < v"1.5"
        cprob = (0,cumsum(collect(fractions))...)
    else
        cprob = (0,cumsum(fractions)...)
    end

    if !(cprob[end] ≈ 1)
        error("fractions do not sum to 1")
    end

    nobs = length(x[1])

    groupindex = zeros(Int,nobs)
    positions = collect(zip(x...))
    c = counter(positions)

    unique_positions = Random.shuffle(rng,collect(keys(c)))

    npos = length(c)
    indices = round.(Int,cprob .* length(c))

    for i=1:length(fractions)
        # set of all position in group i
        set = Set(unique_positions[(indices[i]+1):indices[i+1]])

        for k = 1:nobs
            if positions[k] in set
                groupindex[k] = i
            end
        end
    end
    return groupindex
end
