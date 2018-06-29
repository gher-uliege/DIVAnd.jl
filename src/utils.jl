"""
    cfilled = ufill(c,valex)

Replace values in `c` equal to `valex` by averages of surrounding points.

"""

function ufill(c::Array{T,3},valex::Number) where T
    imax,jmax,kmax = size(c)
    work = zeros(eltype(c),imax+2, jmax+2, kmax+2)
    work2 = zeros(eltype(c),imax+2, jmax+2, kmax+2)

    iwork = zeros(Int8,imax+2, jmax+2, kmax+2)
    iwork2 = zeros(Int8,imax+2, jmax+2, kmax+2)

    cfilled = copy(c)
    ufill!(cfilled,valex,work,work2,iwork,iwork2)

    return cfilled
end



function ufill(c::Array{T,2},valex::Number) where T
    return ufill(reshape(c,(size(c,1), size(c,2), 1)),valex)[:,:,1]
end

"""
    ufill(c::Array{T,2},mask::AbstractArray{Bool}) where T

`mask` is true where `c` is valid.
"""
function ufill(c::Array{T,N},mask::AbstractArray{Bool}) where N where T
    c2 = copy(c)
    # better way
    valex = T(-99999.)
    c2[.!mask] = valex
    
    return ufill(c2,valex)
end

@static if isdefined(:DataArrays)
    ufill(c::DataArray) = ufill(c.data,.!ismissing.(c))
end

function ufill!(c,valexc,work,work2,iwork::Array{Int8,3},iwork2::Array{Int8,3})
    const A1 = 5
    const A2 = 0
    const A3 = 0

    imax,jmax,kmax = size(c)

    for j = 1:jmax+2
        for i = 1:imax+2
            work[i,j,1] = valexc
            iwork[i,j,1] = 0
            work[i,j,kmax+2] = valexc
            iwork[i,j,kmax+2] = 0
        end
    end

    for k = 1:kmax+2
        for i = 1:imax+2
            work[i,1,k] = valexc
            iwork[i,1,k] = 0
            work[i,jmax+2,k] = valexc
            iwork[i,jmax+2,k] = 0
        end
    end

    for k = 1:kmax+2
        for j = 1:jmax+2
            work[1,j,k] = valexc
            iwork[1,j,k] = 0
            work[imax+2,j,k] = valexc
            iwork[imax+2,j,k] = 0
        end
    end

    #
    # copy interior field
    for k = 1:kmax
        for j = 1:jmax
            for i = 1:imax
                work[i+1,j+1,k+1] = c[i,j,k]
                iwork[i+1,j+1,k+1] = 1
                if work[i+1,j+1,k+1] == valexc
                    iwork[i+1,j+1,k+1] = 0
                end
            end
        end
    end

    icount  =  1
    
    while icount > 0
        icount = 0

        for k = 2:kmax+1
            for j = 2:jmax+1
                for i = 2:imax+1

                    work2[i,j,k] = work[i,j,k]
                    iwork2[i,j,k] = iwork[i,j,k]

                    if iwork[i,j,k] == 0
                        work2[i,j,k] = valexc
                        icount = icount+1
                        isom = 0
                        
                        if A1 != 0
                            isom += A1 * (
                                +iwork[i+1,j,k]+iwork[i-1,j,k]
                                +iwork[i,j+1,k]+iwork[i,j-1,k])
                        end
                            
                        if A2 != 0
                            isom += A2 * (
                                iwork[i+1,j+1,k+1]+iwork[i+1,j+1,k-1]
                                +iwork[i+1,j-1,k+1]+iwork[i+1,j-1,k-1]
                                +iwork[i-1,j+1,k+1]+iwork[i-1,j+1,k-1]
                                +iwork[i-1,j-1,k+1]+iwork[i-1,j-1,k-1])
                        end
                        
                        if A3 != 0
                            isom += A3 * (
                                iwork[i,j+1,k+1]+iwork[i,j+1,k-1]
                                + iwork[i,j-1,k+1]+iwork[i,j-1,k-1]
                                + iwork[i+1,j,k+1]+iwork[i+1,j,k-1]
                                + iwork[i-1,j,k+1]+iwork[i-1,j,k-1]
                                + iwork[i+1,j+1,k]+iwork[i+1,j-1,k]
                                + iwork[i-1,j+1,k]+iwork[i-1,j-1,k])
                        end

                        if isom != 0
                            rsom = zero(eltype(c))
                            
                            # interpolate

                            if A1 != 0
                                rsom += A1 * (
                                    +iwork[i+1,j,k]*work[i+1,j,k]
                                    +iwork[i-1,j,k]*work[i-1,j,k]
                                    +iwork[i,j+1,k]*work[i,j+1,k]
                                    +iwork[i,j-1,k]*work[i,j-1,k])
                            end

                            if A2 != 0
                                rsom += A2 * (
                                    iwork[i+1,j+1,k+1]*work[i+1,j+1,k+1]
                                    +iwork[i+1,j+1,k-1]*work[i+1,j+1,k-1]
                                    +iwork[i+1,j-1,k+1]*work[i+1,j-1,k+1]
                                    +iwork[i+1,j-1,k-1]*work[i+1,j-1,k-1]
                                    +iwork[i-1,j+1,k+1]*work[i-1,j+1,k+1]
                                    +iwork[i-1,j+1,k-1]*work[i-1,j+1,k-1]
                                    +iwork[i-1,j-1,k+1]*work[i-1,j-1,k+1]
                                    +iwork[i-1,j-1,k-1]*work[i-1,j-1,k-1])
                            end

                            if A3 != 0                                
                                rsom += A3 * (
                                    iwork[i,j+1,k+1]*work[i,j+1,k+1]
                                    +iwork[i,j+1,k-1]*work[i,j+1,k-1]
                                    +iwork[i,j-1,k+1]*work[i,j-1,k+1]
                                    +iwork[i,j-1,k-1]*work[i,j-1,k-1]
                                    +iwork[i+1,j,k+1]*work[i+1,j,k+1]
                                    +iwork[i+1,j,k-1]*work[i+1,j,k-1]
                                    +iwork[i-1,j,k+1]*work[i-1,j,k+1]
                                    +iwork[i-1,j,k-1]*work[i-1,j,k-1]
                                    +iwork[i+1,j+1,k]*work[i+1,j+1,k]
                                    +iwork[i+1,j-1,k]*work[i+1,j-1,k]
                                    +iwork[i-1,j+1,k]*work[i-1,j+1,k]
                                    +iwork[i-1,j-1,k]*work[i-1,j-1,k])
                            end
                            
                            work2[i,j,k] = rsom/isom
                            iwork2[i,j,k] = 1
                        end
                    end
                end
            end
        end

        for k = 2:kmax+1
            for j = 2:jmax+1
                for i = 2:imax+1
                    work[i,j,k] = work2[i,j,k]
                    iwork[i,j,k] = iwork2[i,j,k]
                end
            end
        end
        #@show icount

    end

# copy interior points
for k = 1:kmax
    for j = 1:jmax
        for i = 1:imax
            c[i,j,k] = work[i+1,j+1,k+1]
        end
    end
end
end




"""
    hx,hy = cgradient(pmn,h)


"""
function cgradient(pmn,h)

    @assert ndims(h) == 2
    
    hx = similar(h)
    hy = similar(h)

    sz = size(h)
    # loop over the domain
    for j = 1:size(h,2)
        # previous j0 and next j1 (but still a valid index)
        j0 = max(j-1,1)
        j1 = min(j+1,sz[2])
        
        for i = 1:size(h,1)
            # previous i0 and next i1 (but still a valid index)
            i0 = max(i-1,1)
            i1 = min(i+1,sz[1])

            # finite difference
            hx[i,j] = (h[i1,j] - h[i0,j]) * pmn[1][i,j]
            hy[i,j] = (h[i,j1] - h[i,j0]) * pmn[2][i,j]

            # centered difference
            if i1 == i0+2
                hx[i,j] = hx[i,j]/2
            end
            
            if j1 == j0+2
                hy[i,j] = hy[i,j]/2
            end

        end
    end
    
    return hx,hy
end


function beforenext(ind,sz::NTuple{N,Int},dim) where N
    # previous and next index (but still a valid index)
    
    ind0 = ntuple(j -> (j == dim ? max(ind[j]-1,1) : ind[j]), N) :: NTuple{N,Int}
    ind1 = ntuple(j -> (j == dim ? min(ind[j]+1,sz[j]) : ind[j]), N) :: NTuple{N,Int}
    return ind0,ind1
end

function cgradient(pmn,h::Array{T,N}, dim) where N where T
    hx = similar(h)
    hy = similar(h)

    sz = size(h)
    
    # loop over the domain
    for ind in CartesianRange(sz)
        #ind = copy(ind)::CartesianIndex{N}
        #ind = ind::CartesianIndex{N}
        
        # previous and next index (but still a valid index)
        ind0,ind1 = beforenext(ind,sz,dim)
        
        #ind0 = ntuple(j -> (j == dim ? max(ind[j]-1,1) : ind[j]), N) :: NTuple{N,Int}
        #ind1 = ntuple(j -> (j == dim ? min(ind[j]+1,sz[j]) : ind[j]), N) :: NTuple{N,Int}        
        
        # finite difference
        hx[ind] = (h[ind1...] - h[ind0...]) * pmn[dim][ind]

        # centered difference
        if (ind[dim] != 1) && (ind[dim] != sz[dim])
            hx[ind] = hx[ind]/2
        end
    end
    
    return hx
end

function cgradientn(pmn,h::Array{T,N}) where N where T
    return ntuple(j -> cgradient(pmn,h,j),N)
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

function lengraddepth(pmn,h::Array{T,2}, L;
                      h2 = h,
                      hmin = 0.001 #m
                      ) where T

    
    # gradient of h
    hx,hy = cgradient(pmn,h)

    normgrad = sqrt.(hx.^2 + hy.^2)

    # avoid divisions by zero
    h2 = max.(h2,hmin)

    # creating the RL field
    RL = 1 ./ (1 + L * normgrad ./ h2)

    #RL[isnan(h)] = valex
    #RL = fill(RL,valex)
    
    return RL
end

function smoothfilter!(x,f::Vector{T},scale) where T
    sz = size(f)
    imax = sz[1]
    
    #      x[1]
    # |     *    |   *   |
    # xf[1]    xf[2]
    
    xf = zeros(T,sz[1]+1)
    xf[1] = x[1] - (x[2]-x[1])/2
    xf[imax+1] = x[imax] + (x[imax]-x[imax-1])/2
        
    for i = 2:imax
        xf[i] = (x[i]+x[i-1])/2
    end

    #@show xf
    ν = sqrt(scale)
    Δx_min = minimum(xf[2:end]-xf[1:end-1])
    #@show Δx_min

    Δt_max = Δx_min^2 / (2 * ν)
    Δt = 0.5 * Δt_max
    
    # 4 t * ν  = 2 scale^2
    # 4 niter * Δt * ν  = 2 scale^2
    # niter = scale^2 / (2 * Δt * ν)
    niter = ceil(Int,scale^2 / (2 * Δt * ν))

    flux = zeros(T,sz[1]+1)
    
    for i = 1:niter
        # flux
        for i = 2:sz[1]
            flux[i] = ν * (f[i] - f[i-1])/(x[i] - x[i-1])
        end

        for i = 1:sz[1]
            f[i] = f[i] + Δt * (flux[i+1] - flux[i])/(xf[i+1] - xf[i])
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

function smoothfilter(x,f::Vector{T},scale) where T
    ff = copy(f)
    smoothfilter!(x,ff,scale)
    return ff
end


"""
    field = divand.random(mask,pmn,len,Nens)

Create `Nens` random fields with the correlation length `len` in 
a domain with the mask `mask` and the metric `pmn`.

See `divand.divandrun` for more information about these parameters.
"""

function random(mask,pmn::NTuple{N,Array{T,N}},len,Nens;
                alpha::Vector{T} = T[],
                moddim::Vector{T} = T[],
                scale_len::Bool = true,
                btrunc = [],
                ) where {N,T}
    
    s = divand.divand_background(
        Val{:sparse},mask,pmn,len,alpha,moddim,scale_len,[];
        btrunc = btrunc);
    
    n = size(s.iB,1)::Int
    z = randn(n,Nens);
    F = cholfact(s.iB::SparseMatrixCSC{T,Int})

    # P pivoting matrix
    # s.iB == P'*L*L'*P
    # F[:UP] ==  L'*P
    
    ff = (F[:UP]) \ z;
    field = divand.unpackens(s.sv,ff)[1] :: Array{T,N+1}
    return field
end


"""
    interp!(xi,fi,x,f)

Interpolate field `fi` (n-dimensional array) defined at `xi` (tuble of
n-dimensional arrays or vectors) onto grid `x` (tuble of n-dimensional arrays).
The interpolated field is stored in `f`.
The grid in `xi` must be align with the axis (e.g. produced by divand.ndgrid).
"""



function interp!(xi::NTuple{N,Vector{T}},
                 fi::Array{T,N},
                 x::NTuple{N,Array{T,Nf}},
                 f::Array{T,Nf}) where {T,N,Nf}   
    itp = interpolate(xi,fi,Gridded(Linear()))

    xpos = zeros(N)
    for i in eachindex(f)
        # position of the i-th location in f
        for j = 1:N
            xpos[j] = x[j][i]
        end
        f[i] = itp[xpos...]
    end
end


function interp!(xi::NTuple{N,Array{T,N}},
                 fi::Array{T,N},
                 x::NTuple{N,Array{T,Nf}},
                 f::Array{T,Nf}) where {T,N,Nf}
    
    # check size
    @assert all([size(xc) == size(fi) for xc in xi])

    # tuple of vector with the varying parts
    xivector = ntuple(j -> xi[j][[(i==j ? (:) : 1 ) for i in 1:N]...], N) :: NTuple{N,Vector{T}}
    interp!(xivector,fi,x,f)
end

"""
    f = interp(xi,fi,x)

Interpolate field `fi` (n-dimensional array) defined at `xi` (tuble of
n-dimensional arrays or vectors) onto grid `x` (tuble of n-dimensional arrays).
The grid in `xi` must be align with the axis (e.g. produced by divand.ndgrid).
"""

function interp(xi,fi,x)
    f = similar(x[1])
    interp!(xi,fi,x,f)
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

function backgroundfile(fname,varname)
    ds = Dataset(fname)
    lon = nomissing(ds["lon"][:])
    lat = nomissing(ds["lat"][:])
    depth = nomissing(ds["depth"][:])
    #time = ds["time"][:].data
    
    v = ds[varname]
    x = (lon,lat,depth)    
    
    return function (xi,n,value,trans)
        
        vn = zeros(size(v[:,:,:,n]))
        vn .= map((x -> ismissing(x) ? NaN : x), v[:,:,:,n]);

        
        vn .= trans.(divand.ufill(vn,.!isnan.(vn)))
        fi = divand.interp(x,vn,xi)

        return vn,value - fi
    end
end


dayssince(dt; t0 = DateTime(1900,1,1)) = Dates.value.(dt - t0)/1000/60/60/24;
