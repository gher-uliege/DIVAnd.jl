
"""
    checkresolution(mask,pmn,len)

Returns a warning of the resolution is too coarse relative to the correlation
length. The resolution must be at least 2 times finer than the correlation
length.
"""
function checkresolution(
    mask,
    pmn::NTuple{N,Array{T1,N}},
    len::NTuple{N,Array{T2,N}},
) where {N,T1,T2}
    for i = 1:length(pmn)
        for j in CartesianIndices(pmn[i])
            if ((pmn[i][j] * len[i][j] <= 2) && mask[j]) && (len[i][j] != 0.0)
                res = 1 / pmn[i][j]
                @warn "resolution ($res) is too coarse for correlation length $(len[i][j]) in dimension $i at indices $j (skipping further tests). It is recommended that the resolution is at least 2 times finer than the correlation length."
                #error("stop")
                break
            end
        end
    end
end

checkresolution(mask, pmn, len) = checkresolution(mask, pmn, len_harmonize(len, mask))


function checkdepth(depthr)
    if length(unique(depthr)) !== length(depthr)
        error("Depth levels should be unique $(depth)")
    end

    if any(depthr[2:end] .<= depthr[1:end-1])
        error("Depth levels should increase monotonically")
    end
end

"""
    cfilled = ufill(c,valex)

Replace values in `c` equal to `valex` by averages of surrounding points.
`valex` should not be NaN; use `ufill(c,isfinite.(c))` or
`ufill(c,.!isnan.(c))` instead.
"""
function ufill(c::Array{T,3}, valex::Number) where {T}
    #JLD2.@save "/tmp/tmp-ufill.jld2" c valex

    imax, jmax, kmax = size(c)
    work = zeros(eltype(c), imax + 2, jmax + 2, kmax + 2)
    work2 = zeros(eltype(c), imax + 2, jmax + 2, kmax + 2)

    iwork = zeros(Int8, imax + 2, jmax + 2, kmax + 2)
    iwork2 = zeros(Int8, imax + 2, jmax + 2, kmax + 2)

    cfilled = copy(c)

    if any(sum(c .!= valex, dims = [1, 2]) .== 0)
        minval, minindex = findmin(sum(c .!= valex, dims = [1, 2])[:])
        @show valex
        @show size(c)
        error("some slices completely masked: k = $(minindex) of array $(size(c))")
    end
    ufill!(cfilled, valex, work, work2, iwork, iwork2)

    return cfilled
end



function ufill(c::Array{T,2}, valex::Number) where {T}
    return ufill(reshape(c, (size(c, 1), size(c, 2), 1)), valex)[:, :, 1]
end

function ufill(c::Array{T,1}, valex::Number) where {T}
    return ufill(reshape(c, (:, 1, 1)), valex)[:, 1, 1]
end

"""
    ufill(c::Array{T,2},mask::AbstractArray{Bool}) where T

`mask` is true where `c` is valid.
"""
function ufill(c::Array{T,N}, mask::AbstractArray{Bool}) where {N} where {T}
    c2 = copy(c)
    # better way
    valex = T(-99999.)
    c2[.!mask] .= valex

    return ufill(c2, valex)
end

function ufill!(c, valexc, work, work2, iwork::Array{Int8,3}, iwork2::Array{Int8,3})
    A1 = 5
    A2 = 0
    A3 = 0

    imax, jmax, kmax = size(c)

    for j = 1:jmax+2
        for i = 1:imax+2
            work[i, j, 1] = valexc
            iwork[i, j, 1] = 0
            work[i, j, kmax+2] = valexc
            iwork[i, j, kmax+2] = 0
        end
    end

    for k = 1:kmax+2
        for i = 1:imax+2
            work[i, 1, k] = valexc
            iwork[i, 1, k] = 0
            work[i, jmax+2, k] = valexc
            iwork[i, jmax+2, k] = 0
        end
    end

    for k = 1:kmax+2
        for j = 1:jmax+2
            work[1, j, k] = valexc
            iwork[1, j, k] = 0
            work[imax+2, j, k] = valexc
            iwork[imax+2, j, k] = 0
        end
    end

    #
    # copy interior field
    for k = 1:kmax
        for j = 1:jmax
            for i = 1:imax
                work[i+1, j+1, k+1] = c[i, j, k]
                iwork[i+1, j+1, k+1] = 1
                if work[i+1, j+1, k+1] == valexc
                    iwork[i+1, j+1, k+1] = 0
                end
            end
        end
    end

    icount = 1

    while icount > 0
        icount = 0

        for k = 2:kmax+1
            for j = 2:jmax+1
                for i = 2:imax+1

                    work2[i, j, k] = work[i, j, k]
                    iwork2[i, j, k] = iwork[i, j, k]

                    if iwork[i, j, k] == 0
                        work2[i, j, k] = valexc
                        icount = icount + 1
                        isom = 0

                        if A1 != 0
                            isom += A1 *
                                    (+iwork[i+1, j, k] + iwork[i-1, j, k] +
                                     iwork[i, j+1, k] + iwork[i, j-1, k])
                        end

                        if A2 != 0
                            isom += A2 *
                                    (iwork[i+1, j+1, k+1] + iwork[i+1, j+1, k-1] +
                                     iwork[i+1, j-1, k+1] + iwork[i+1, j-1, k-1] +
                                     iwork[i-1, j+1, k+1] + iwork[i-1, j+1, k-1] +
                                     iwork[i-1, j-1, k+1] + iwork[i-1, j-1, k-1])
                        end

                        if A3 != 0
                            isom += A3 *
                                    (iwork[i, j+1, k+1] + iwork[i, j+1, k-1] +
                                     iwork[i, j-1, k+1] + iwork[i, j-1, k-1] +
                                     iwork[i+1, j, k+1] + iwork[i+1, j, k-1] +
                                     iwork[i-1, j, k+1] + iwork[i-1, j, k-1] +
                                     iwork[i+1, j+1, k] + iwork[i+1, j-1, k] +
                                     iwork[i-1, j+1, k] + iwork[i-1, j-1, k])
                        end

                        if isom != 0
                            rsom = zero(eltype(c))

                            # interpolate

                            if A1 != 0
                                rsom += A1 *
                                        (+iwork[i+1, j, k] * work[i+1, j, k] +
                                         iwork[i-1, j, k] * work[i-1, j, k] +
                                         iwork[i, j+1, k] * work[i, j+1, k] +
                                         iwork[i, j-1, k] * work[i, j-1, k])
                            end

                            if A2 != 0
                                rsom += A2 *
                                        (iwork[i+1, j+1, k+1] * work[i+1, j+1, k+1] +
                                         iwork[i+1, j+1, k-1] * work[i+1, j+1, k-1] +
                                         iwork[i+1, j-1, k+1] * work[i+1, j-1, k+1] +
                                         iwork[i+1, j-1, k-1] * work[i+1, j-1, k-1] +
                                         iwork[i-1, j+1, k+1] * work[i-1, j+1, k+1] +
                                         iwork[i-1, j+1, k-1] * work[i-1, j+1, k-1] +
                                         iwork[i-1, j-1, k+1] * work[i-1, j-1, k+1] +
                                         iwork[i-1, j-1, k-1] * work[i-1, j-1, k-1])
                            end

                            if A3 != 0
                                rsom += A3 *
                                        (iwork[i, j+1, k+1] * work[i, j+1, k+1] +
                                         iwork[i, j+1, k-1] * work[i, j+1, k-1] +
                                         iwork[i, j-1, k+1] * work[i, j-1, k+1] +
                                         iwork[i, j-1, k-1] * work[i, j-1, k-1] +
                                         iwork[i+1, j, k+1] * work[i+1, j, k+1] +
                                         iwork[i+1, j, k-1] * work[i+1, j, k-1] +
                                         iwork[i-1, j, k+1] * work[i-1, j, k+1] +
                                         iwork[i-1, j, k-1] * work[i-1, j, k-1] +
                                         iwork[i+1, j+1, k] * work[i+1, j+1, k] +
                                         iwork[i+1, j-1, k] * work[i+1, j-1, k] +
                                         iwork[i-1, j+1, k] * work[i-1, j+1, k] +
                                         iwork[i-1, j-1, k] * work[i-1, j-1, k])
                            end

                            work2[i, j, k] = rsom / isom
                            iwork2[i, j, k] = 1
                        end
                    end
                end
            end
        end

        for k = 2:kmax+1
            for j = 2:jmax+1
                for i = 2:imax+1
                    work[i, j, k] = work2[i, j, k]
                    iwork[i, j, k] = iwork2[i, j, k]
                end
            end
        end
        #@show icount

    end

# copy interior points
    for k = 1:kmax
        for j = 1:jmax
            for i = 1:imax
                c[i, j, k] = work[i+1, j+1, k+1]
            end
        end
    end
end

"""
    directions = vonNeumannNeighborhood(mask)

Return a vector will all search directions corresponding to the Von Neumann
neighborhood in N dimensions where N is the dimension of the boolean array
`mask`.
"""
function vonNeumannNeighborhood(mask::AbstractArray{Bool,N}) where {N}
    return [CartesianIndex(ntuple(i -> (i == j ? s : 0), Val(N))) for j = 1:N for s in [
        -1,
        1,
    ]]
end

"""
    m = floodfillpoint(mask,I,directions)

Fill the binary mask starting at index `I` (`CartesianIndex`). All element
directly connected to the starting location `I` will be `true` without crossing
any element equal to `false` in `mask`. Per default the value of `I` is the
first true element in `mask` and `directions ` correspond to the Von Neumann
neighborhood.
"""
function floodfillpoint(
    mask,
    I = findfirst(mask),
    directions = vonNeumannNeighborhood(mask),
)
    m = falses(size(mask))

    m[I] = true

    anyflip = true

    CI = CartesianIndices(size(m))

    while anyflip
        anyflip = false

        for I in CI
            if m[I]
                for dir in directions

                    i1 = I + dir
                    if checkbounds(Bool, m, i1)
                        if mask[i1] && !m[i1]
                            m[i1] = true
                            anyflip = true
                        end
                    end
                end
            end
        end

        if !anyflip
            break
        end
    end
    return m
end

"""
    label = floodfill(mask)

Attribute an integer number (a numeric label) to every element in mask such that
all grid points connected by a von Neumann neighborhood (without crossing
elements which are `false` in mask) have the same label.
Labels are sorted such that the label 1 corresponds to the largest area, label 2
the 2nd largest and so on.
"""
function floodfill(mask, directions = vonNeumannNeighborhood(mask))
    m = copy(mask)
    index = zeros(Int, size(mask))
    area = Int[]

    l = 0
    while any(m)
        l = l + 1
        ml = floodfillpoint(m, findfirst(m), directions)
        index[ml] .= l
        m[ml] .= false
        push!(area, sum(ml))
    end

    sortp = sortperm(area; rev = true)
    for I in eachindex(index)
        tmp = index[I]
        if tmp != 0
            index[I] = sortp[tmp]
        end
    end

    return index
end

"""
    hx,hy = cgradient(pmn,h)

"""
function cgradient(pmn, h)

    @assert ndims(h) == 2

    hx = similar(h)
    hy = similar(h)

    sz = size(h)
    # loop over the domain
    for j = 1:size(h, 2)
        # previous j0 and next j1 (but still a valid index)
        j0 = max(j - 1, 1)
        j1 = min(j + 1, sz[2])

        for i = 1:size(h, 1)
            # previous i0 and next i1 (but still a valid index)
            i0 = max(i - 1, 1)
            i1 = min(i + 1, sz[1])

            # finite difference
            hx[i, j] = (h[i1, j] - h[i0, j]) * pmn[1][i, j]
            hy[i, j] = (h[i, j1] - h[i, j0]) * pmn[2][i, j]

            # centered difference
            if i1 == i0 + 2
                hx[i, j] = hx[i, j] / 2
            end

            if j1 == j0 + 2
                hy[i, j] = hy[i, j] / 2
            end

        end
    end

    return hx, hy
end


function beforenext(ind, sz::NTuple{N,Int}, dim) where {N}
    # previous and next index (but still a valid index)

    ind0 = ntuple(j -> (j == dim ? max(ind[j] - 1, 1) : ind[j]), N)::NTuple{N,Int}
    ind1 = ntuple(j -> (j == dim ? min(ind[j] + 1, sz[j]) : ind[j]), N)::NTuple{N,Int}
    return ind0, ind1
end

function cgradient(pmn, h::Array{T,N}, dim) where {N} where {T}
    hx = similar(h)
    hy = similar(h)

    sz = size(h)

    cranges = CartesianIndices(ntuple(i -> 1:size(h, i), N))

    # loop over the domain
    for ind in cranges
        #ind = copy(ind)::CartesianIndex{N}
        #ind = ind::CartesianIndex{N}

        # previous and next index (but still a valid index)
        ind0, ind1 = beforenext(ind, sz, dim)

        #ind0 = ntuple(j -> (j == dim ? max(ind[j]-1,1) : ind[j]), N) :: NTuple{N,Int}
        #ind1 = ntuple(j -> (j == dim ? min(ind[j]+1,sz[j]) : ind[j]), N) :: NTuple{N,Int}

        # finite difference
        hx[ind] = (h[ind1...] - h[ind0...]) * pmn[dim][ind]

        # centered difference
        if (ind[dim] != 1) && (ind[dim] != sz[dim])
            hx[ind] = hx[ind] / 2
        end
    end

    return hx
end

function cgradientn(pmn, h::Array{T,N}) where {N} where {T}
    return ntuple(j -> cgradient(pmn, h, j), N)
end


"""
    RL = lengraddepth(pmn,h, L;
                      h2 = h,
                      hmin = 0.001
                      )

Create the relative correlation length-scale field `RL` based on the bathymetry
`h` and the metric `pmn` (tuple of arrays). Effectively the correlation-length
scale is close to zero if the relative bathymetry gradients (|∇h|/h) are smaller
 than the length-scale `L` (in consistent units as `pmn`).

R_L = 1 / (1 + L |∇h| / max(h2,hmin))

Per default `h2` is equal to `h`. The depth `h` must be positive. `hmin` must
have the same units as h (usually meters).
"""
function lengraddepth(
    pmn,
    h::Array{T,2},
    L;
    h2 = h,
    hmin = 0.001, #m
) where {T}


    # gradient of h
    hx, hy = cgradient(pmn, h)

    normgrad = sqrt.(hx.^2 + hy.^2)

    # avoid divisions by zero
    h2 = max.(h2, hmin)

    # creating the RL field
    RL = 1 ./ (1 .+ L * normgrad ./ h2)

    #RL[isnan(h)] = valex
    #RL = fill(RL,valex)

    return RL
end

function smoothfilter!(x, f::Vector{T}, scale) where {T}
    sz = size(f)
    imax = sz[1]

    #      x[1]
    # |     *    |   *   |
    # xf[1]    xf[2]

    xf = zeros(T, sz[1] + 1)
    xf[1] = x[1] - (x[2] - x[1]) / 2
    xf[imax+1] = x[imax] + (x[imax] - x[imax-1]) / 2

    for i = 2:imax
        xf[i] = (x[i] + x[i-1]) / 2
    end

    #@show xf
    ν = sqrt(scale)
    Δx_min = minimum(xf[2:end] - xf[1:end-1])
    #@show Δx_min

    Δt_max = Δx_min^2 / (2 * ν)
    Δt = 0.5 * Δt_max

    # 4 t * ν  = 2 scale^2
    # 4 niter * Δt * ν  = 2 scale^2
    # niter = scale^2 / (2 * Δt * ν)
    niter = ceil(Int, scale^2 / (2 * Δt * ν))

    flux = zeros(T, sz[1] + 1)

    for i = 1:niter
        # flux
        for i = 2:sz[1]
            flux[i] = ν * (f[i] - f[i-1]) / (x[i] - x[i-1])
        end

        for i = 1:sz[1]
            f[i] = f[i] + Δt * (flux[i+1] - flux[i]) / (xf[i+1] - xf[i])
        end
    end

end


"""
    ff = smoothfilter(x,f,scale)

Smooth the function `f` defined on `x` by solving the diffusion equation

∂ₜ ϕ = ν ∂²ₓ ϕ

`scale` is the spatial scales of the removed length-scales.
It is defined as 2Tν  where T is the integration time.

It uses the Greens functions for 1D diffusion:
1/sqrt(4 π ν t) * exp(-x^2 / (4νt))

"""
function smoothfilter(x, f::Vector{T}, scale) where {T}
    ff = copy(f)
    smoothfilter!(x, ff, scale)
    return ff
end

"""
All weights have to be positive (and different from zero).
"""
function smoothfilter_weighted(x, f::Vector{T}, w, scale) where {T}
    wff = smoothfilter(x,w.*f,scale)
    wf = smoothfilter(x,w,scale)
    return wff ./ wf, wf
end


"""
    field = DIVAnd.random(mask,pmn,len,Nens)

Create `Nens` random fields with the correlation length `len` in
a domain with the mask `mask` and the metric `pmn`.

See `DIVAnd.DIVAndrun` for more information about these parameters.
"""
function random(
    mask,
    pmn::NTuple{N,Array{T,N}},
    len,
    Nens;
    alpha::Vector{T} = T[],
    moddim::Vector{T} = T[],
    scale_len::Bool = true,
    btrunc = [],
) where {N,T}

    s = DIVAnd.DIVAnd_background(
        Val{:sparse},
        mask,
        pmn,
        len,
        alpha,
        moddim,
        scale_len,
        [];
        btrunc = btrunc,
    )

    n = size(s.iB, 1)::Int
    z = randn(n, Nens)

    F = cholesky(s.iB::SparseMatrixCSC{T,Int})
    F_UP = F.UP

    # P pivoting matrix
    # s.iB == P'*L*L'*P
    # F[:UP] ==  L'*P

    ff = F_UP \ z
    field = DIVAnd.unpackens(s.sv, ff)[1]::Array{T,N + 1}
    return field
end


"""
    interp!(xi,fi,x,f)

Interpolate field `fi` (n-dimensional array) defined at `xi` (tuble of
n-dimensional arrays or vectors) onto grid `x` (tuble of n-dimensional arrays).
The interpolated field is stored in `f`.
The grid in `xi` must be align with the axis (e.g. produced by DIVAnd.ndgrid).
"""
function interp!(
    xi::NTuple{N,Vector{T}},
    fi::Array{T,N},
    x::NTuple{N,Array{T,Nf}},
    f::Array{T,Nf},
) where {T,N,Nf}

    # https://github.com/JuliaMath/Interpolations.jl/issues/237
    itp = extrapolate(interpolate(xi, fi, Gridded(Linear())), Line())

    xpos = zeros(N)
    for i in eachindex(f)
        # position of the i-th location in f
        for j = 1:N
            xpos[j] = x[j][i]
        end
        f[i] = itp(xpos...)
    end
end


function interp!(
    xi::NTuple{N,Array{T,N}},
    fi::Array{T,N},
    x::NTuple{N,Array{T,Nf}},
    f::Array{T,Nf},
) where {T,N,Nf}

    # check size
    @assert all([size(xc) == size(fi) for xc in xi])

    # tuple of vector with the varying parts
    xivector = ntuple(j -> xi[j][[(i == j ? (:) : 1) for i = 1:N]...], N)::NTuple{
        N,
        Vector{T},
    }
    interp!(xivector, fi, x, f)
end

"""
    f = interp(xi,fi,x)

Interpolate field `fi` (n-dimensional array) defined at `xi` (tuble of
n-dimensional arrays or vectors) onto grid `x` (tuble of n-dimensional arrays).
The grid in `xi` must be align with the axis (e.g. produced by DIVAnd.ndgrid).
"""
function interp(xi, fi, x)
    f = similar(x[1])
    interp!(xi, fi, x, f)
    return f
end




"""
    fun = backgroundfile(fname,varname)

Return a function `fun` which is used in DIVAnd to make
anomalies out of observations based relative to the field
defined in the NetCDF variable `varname` in the NetCDF file
`fname`. It is assumed that the NetCDF variables has the variable
`lon`, `lat` and `depth`. And that the NetCDF variable is defined on the
same grid as the analysis.

"""
function backgroundfile(fname, varname)
    ds = Dataset(fname)
    lon = nomissing(ds["lon"][:])
    lat = nomissing(ds["lat"][:])
    depth = nomissing(ds["depth"][:])
    #time = ds["time"][:].data

    v = ds[varname]
    x = (lon, lat, depth)

    return function (xi, n, value, trans; selection = [], obstime = nothing)

        vn = zeros(size(v[:, :, :, n]))
        vn .= map((x -> ismissing(x) ? NaN : x), v[:, :, :, n])


        vn .= trans.(DIVAnd.ufill(vn, .!isnan.(vn)))
        fi = DIVAnd.interp(x, vn, xi)

        return vn, value - fi
    end
end

"""
    fun = backgroundfile(fname,varname,TS)

Return a function `fun` which is used in DIVAnd to make
anomalies out of observations based relative to the field
defined in the NetCDF variable `varname` in the NetCDF file
`fname`. It is assumed that the NetCDF variables has the variable
`lon`, `lat` and `depth`. And that the NetCDF variable is defined on the
same grid as the analysis and was generated according to the provided time selector
`TS` (TimeSelectorYearListMonthList or TimeSelectorRunningAverage).

!!! note

    At all vertical levels, there should at least one sea point.
"""
function backgroundfile(
    fname,
    varname,
    TS::AbstractTimeSelector
)

    ds = Dataset(fname)
    lon = nomissing(ds["lon"][:])
    lat = nomissing(ds["lat"][:])
    depth = nomissing(ds["depth"][:])

    v = ds[varname]
    x = (lon, lat, depth)
    TSbackground = TS

    return function (xi, n, value, trans; selection = [], obstime = nothing)
        # check which background estimate has the best overlap
        overlap = zeros(Int, length(TSbackground))
        for timeindex = 1:length(TSbackground)
            sel = select(TSbackground, timeindex, obstime)
            overlap[timeindex] = sum(selection .& sel)
        end

        nbackground = findmax(overlap)[2]

        @info "analysis time index $n uses the backgrond time index $nbackground"

        vn = zeros(size(v[:, :, :, nbackground]))
        vn .= map((x -> ismissing(x) ? NaN : x), v[:, :, :, nbackground])

        vn .= trans.(DIVAnd.ufill(vn, .!isnan.(vn)))
        fi = DIVAnd.interp(x, vn, xi)

        return vn, value - fi
    end
end

"""
    dayssince(dt; t0 = DateTime(1900,1,1))

Number of days since a starting day `t0` (1900-01-01 per default).
"""
dayssince(dt; t0 = DateTime(1900, 1, 1)) = Dates.value.(dt - t0) / 1000 / 60 / 60 / 24;




function _diffusionfix!(ivol, nus, α, nmax, x0, x)
    work1 = similar(x)
    x[:] = x0

    for niter = 1:nmax
        DIVAnd.DIVAnd_laplacian_apply!(ivol, nus, x, work1)
        for i = 1:length(x0)
            if x0[i] != 0
                x[i] = x[i] + α * work1[i]
            end
        end
    end

end


"""
    mergedfield = hmerge(field,L)

Merge several `field[:,:,1]`, `field[:,:,2]`,... into a single 2d field
`mergedfield` values equal to NaN are ignored. This function is typically used
to merge different DIVAnd anayses.
"""
function hmerge(f, L)
    # L ∼ (α * nmax)²
    # nmax ∼ √(L)/α

    weight0 = Float64.(isfinite.(f))

    mask, pmn = DIVAnd.DIVAnd_rectdom(1:size(f, 1), 1:size(f, 2))
    ivol, nus = DIVAnd.DIVAnd_laplacian_prepare(
        mask,
        pmn,
        (ones(size(mask)), ones(size(mask))),
    )

    α = 0.1
    nmax = round(Int, sqrt(L) / α)
    @debug "nmax: $(nmax)"

    #nmax = 20;
    weight = similar(weight0)

    for k = 1:size(weight, 3)
        wk0 = @view weight0[:, :, k]
        wk = @view weight[:, :, k]

        _diffusionfix!(ivol, nus, α, nmax, wk0, wk)
    end
    f[.!isfinite.(f)] .= 0

    weight = weight.^2

    f2 = (sum(weight .* f, dims = 3)./sum(weight, dims = 3))[:, :, 1]

    return f2
end



"""
    fun = velocityfile(fname,(varnameu,varnamev),TSvelocity,scale)
    fun = velocityfile(fname,(varnameu,varnamev,varnamew),TSvelocity,scale)

Return a function `fun` which is used in DIVAnd as a advection constraint using
fields defined in the NetCDF variable `varnameu`, `varnamev` and `varnamew` 
(zonal, meridional and vertical velocity components) in the NetCDF file
`fname`. If the parameter `varnamew` is omitted, the vertical velicity is neglected.
It is assumed that the NetCDF variables has the variable
`lon`, `lat` and `depth` and that the fields have been average according to the provided time selector
`TSvelocity` (`TimeSelectorYearListMonthList` or `TimeSelectorRunningAverage`).

See also `DIVAnd.average_files`.
!!! note

    NetCDF _FillValues are treated as zeros.
"""
function velocityfile(
    fname,
    varnames,
    TSvelocity:: DIVAnd.AbstractTimeSelector,
    scale
)

    ds = Dataset(fname)
    lon = nomissing(ds["lon"][:])
    lat = nomissing(ds["lat"][:])
    depth = nomissing(ds["depth"][:])


    return function (xi::NTuple{N,AbstractArray{T,N}}, veltime) where N where T
        # check which velocity estimate has the best overlap
        overlap = zeros(Int, length(TSvelocity))
        for timeindex = 1:length(TSvelocity)
            sel = DIVAnd.select(TSvelocity, timeindex, veltime)
            overlap[timeindex] = sum(sel)
        end

        nvelocity = findmax(overlap)[2]

        @info "analysis uses the velocity time index $nvelocity"

        x = (T.(lon), T.(lat), T.(depth))

        u = ds[varnames[1]]
        tmpu = u[:, :, :, nvelocity]
        un = T.(replace(tmpu, missing => zero(T)))
        ui = DIVAnd.interp(x, un, xi)

        v = ds[varnames[2]]
        tmpv = v[:, :, :, nvelocity]
        vn = T.(replace(tmpv, missing => zero(T)))
        vi = DIVAnd.interp(x, vn, xi)

        wi =
            if length(varnames) == 3
                w = ds[varnames[3]]
                tmpw = w[:, :, :, nvelocity]
                wn = T.(replace(tmpw, missing => zero(T)))
                DIVAnd.interp(x, vn, xi)
            else
                zeros(size(ui))
            end

        return (scale*ui,scale*vi,scale*wi)
    end
end



"""
    meanerr(obsconstrain)

Compute the mean error variance of the observational contraints `obsconstrain`
ignoring values equal to `Inf` (due to e.g. observations outside of the grid).
"""
meanerr(obsconstrain) = mean(diag(obsconstrain.R)[isfinite.(diag(obsconstrain.R))])


"""
    average_files(filenames,varnameu::AbstractString,TSvelocity,outfilename,outvarname::AbstractString)

Averaged the gridded variables `varnameu` from the NetCDF files from `filenames`
over time following the time selector `TSvelocity`.
"""
function average_files(filenames,varnameu::AbstractString,TSvelocity,outfilename,outvarname::AbstractString)
    ds = Dataset(filenames; aggdim = "time")

    lon = ds["lon"][:];
    lat = ds["lat"][:];
    depth = ds["depth"][:];
    time = ds["time"][:];

    u = ds[varnameu];

    @info("average $varnameu")
    meanu = average_field(time,u,TSvelocity)

    close(ds)

    mode = "c"
    if isfile(outfilename)
        mode = "a"
    end

    fillvalue = NCDatasets.NC_FILL_FLOAT

    Dataset(outfilename,mode, attrib = [
        "history"                     => "generated by $(@__FILE__)",
    ]) do ds

        if mode == "c"
            # Dimensions

            ds.dim["time"] = size(meanu,4)
            ds.dim["depth"] = size(meanu,3)
            ds.dim["lat"] = size(meanu,2)
            ds.dim["lon"] = size(meanu,1)

            # Declare variables
            ncdepth = defVar(ds,"depth", Float32, ("depth",), attrib = [
                "units"                     => "m",
                "long_name"                 => "depth",
                "standard_name"             => "depth",
                "positive"                  => "down",
            ])

            nclon = defVar(ds,"lon", Float32, ("lon",), attrib = [
                "units"                     => "degrees_east",
                "long_name"                 => "longitude",
                "standard_name"             => "longitude",
            ])

            nclat = defVar(ds,"lat", Float32, ("lat",), attrib = [
                "units"                     => "degrees_north",
                "long_name"                 => "latitude",
                "standard_name"             => "latitude",
            ])

            # nctime = defVar(ds,"time", Float64, ("time",), attrib = [
            #     "units"                     => "days since 1970-01-01 00:00:00",
            #     "long_name"                 => "time",
            #     "standard_name"             => "time",
            #     "calendar"                  => "standard",
            # ])


            nclon[:] = lon
            nclat[:] = lat
            ncdepth[:] = depth
            #nctime[:] = time
        end

        ncu = defVar(ds,outvarname, Float32, ("lon", "lat", "depth", "time"), attrib = [
            "_FillValue"                => fillvalue,
        ])
        ncu[:,:,:,:] = replace(meanu,NaN => missing)

    end # closing ds
    return nothing
end


function average_files(filenames,varnames::Tuple,TSvelocity,outfilename,outvarnames::Tuple)
    for (varname,outvarname) in zip(varnames,outvarnames)
        average_files(filenames,varname,TSvelocity,outfilename,outvarname)
    end
end

function average_field(time,u,TSvelocity)
    sz = size(u)
    meanu = zeros(sz[1],sz[2],sz[3],length(TSvelocity))

    for n = 1:length(TSvelocity)
        @info("time instance $n of $(length(TSvelocity))")
        count = 0

        for n2 in findall(DIVAnd.select(TSvelocity,n,time))
            meanu[:,:,:,n] += nomissing(u[:,:,:,n2],NaN)
            count += 1
        end

        meanu[:,:,:,n] = meanu[:,:,:,n]/count
    end
    return meanu
end
