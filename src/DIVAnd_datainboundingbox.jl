"""
    xn,fn,indexes,Rn = DIVAnd_datainboundingbox(xi,x,f;Rmatrix=())

# Input:
  `xi`: tuple with n elements. Every element represents a coordinate
  of the final grid on which the observations are interpolated

* `x`: tuple with n elements. Every element represents a coordinate of
  the observations

* `f`: value of the observations

* `Rmatrix`: error variance of the observations (normalized by the error variance of the background field). `epsilon2` can be a scalar (all observations have the same error variance and their errors are decorrelated), a vector (all observations can have a difference error variance and their errors are decorrelated) or a matrix (all observations can have a difference error variance and their errors can be correlated). If `epsilon2` is a scalar, it is thus the *inverse of the signal-to-noise ratio*.



# Output:

`xn`: tuple with n elements. Every element represents a coordinate of
  the observations which falls in the bounding box defined by xi
`fn`: the corresponding data
`indexes:` the indexes in the original array retained
`Rn`: the new error variance

"""
function DIVAnd_datainboundingbox(xi, x, f; Rmatrix = ())

    n = length(xi)

    maxxi = zeros(n)
    minxi = zeros(n)

    sel = trues(size(x[1], 1))

    for i = 1:n
        maxxi[i] = maximum(xi[i])
        minxi[i] = minimum(xi[i])
    end

    for j = 1:size(x[1], 1)
        for i = 1:n
            # JMB: To avoid a point taken by tile be on the boundary and not seen by divandrun to produce a NaN, use a slightly severe test ?
            # Alex: Can be too strict for surface data. Maybe it is not necessary anymore after commit 8af890d1cf06fda8f2520bb01a28b5b5128051bc
            if !(minxi[i] <= x[i][j] <= maxxi[i])
            #dxx = (maxxi[i] - minxi[i]) * 0.0000000001
            #if !(minxi[i] + dxx < x[i][j] < maxxi[i] - dxx)
                sel[j] = false
                break
            end
        end
    end

    if Rmatrix == ()
        Rm = ()
    else
        if isa(Rmatrix, Number)
            Rm = Rmatrix
        else
            if ndims(Rmatrix) == 1
                Rm = Rmatrix[sel]
            else
                Rm = Rmatrix[sel, sel]
            end

        end
    end

    return ([xx[sel] for xx in x]...,), f[sel], findall(sel), Rm
end
