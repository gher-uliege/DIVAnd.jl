


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

    for i = 1:length(uniquepos)
        for j = 1:length(allpos)
            if uniquepos[i] == allpos[j]
                count[i] += 1
                sumval[i] += val[j]
                sumval2[i] += val[j]^2
            end
        end
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
