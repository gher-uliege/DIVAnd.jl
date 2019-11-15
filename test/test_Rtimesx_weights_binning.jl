import DIVAnd
using Test
using DelimitedFiles
using DataStructures
using Missings
using NCDatasets
using Interpolations
using Random
using Statistics

#
varname = "Salinity"
filename = "WOD-Salinity.nc"


bathname = joinpath(
    dirname(@__FILE__),
    "..",
    "..",
    "DIVAnd-example-data",
    "Global",
    "Bathymetry",
    "gebco_30sec_16.nc",
)
bathisglobal = true

obsname = joinpath(
    dirname(@__FILE__),
    "..",
    "..",
    "DIVAnd-example-data",
    "Provencal",
    "WOD-Salinity.nc",
)

cdilist = joinpath(dirname(@__FILE__), "..", "data", "CDI-list-export.csv")


if !isfile(bathname)
    @info("download bathymetry $bathname")
    bathname = download("https://dox.ulg.ac.be/index.php/s/U0pqyXhcQrXjEUX/download")
end


if !isfile(obsname)
    @info("download observations $obsname")
    obsname = download("https://dox.ulg.ac.be/index.php/s/PztJfSEnc8Cr3XN/download")
end

obsvalue, obslon, obslat, obsdepth, obstime, obsids = DIVAnd.loadobs(
    Float64,
    obsname,
    "Salinity",
)



function binning(gridx::Tuple, x, v)
    # unstaggered coordinate
    gridx_unstagger = DIVAnd.unstagger.(gridx)
    sz = length.(gridx)

    vb = zeros(eltype(v), sz)
    count = zeros(Int, sz)
    nout = 0
    n = length(gridx)

    for j = 1:length(v)
        ind = CartesianIndex(ntuple(l -> DIVAnd.findin(gridx_unstagger[l], x[l][j]), Val(n)))

        if checkbounds(Bool,vb,ind)
            vb[ind] = vb[ind] + v[j]
            count[ind] += 1
        else
            nout += 1
        end
    end

    return vb ./ count, count, vb, nout
end


function weight_RtimesOne_binning(x, len)
    n = length(x)

    # grid finer than a factor of 10
    dx = ntuple(i -> len[i]/10,Val(n))

    gridx = ntuple(i -> minimum(x[i]):dx[i]:maximum(x[i])+dx[i], Val(n))
    sz = length.(gridx)

    v = ones(size(x[1]))
    m2, count2, vb2, nout2 = binning(gridx, x, v)

    mask = trues(sz)
    pmn = ntuple(i -> ones(sz)/dx[1], Val(n))

    c0 = Float64.(count2);
    c = zeros(size(c0))

    # adusted length
    coef = sqrt(2)
    len_adjusted = ntuple(i -> coef * len[i], Val(n))

    # "diffusion" coefficient
    nu = ntuple(i -> fill(len_adjusted[i].^2,sz),Val(n))

    # compute inverses of cell volumne and staggered scaled coefficients
    ivol, nus = DIVAnd.DIVAnd_laplacian_prepare(mask,pmn,nu)

    # maximum allowed time step
    α0 = 1 / (2 * sum(ntuple(i -> maximum(pmn[i].^2 .* nu[i]),Val(n))))

    # 10% safety margin
    α = α0 / 1.1

    # number of iterations 1/(2*α) (rounded)
    nmax = round(Int, 1 / (2 * α))

    # 4* L² α*nmax ≈ 2 L² = L'²
    @debug "α0: $α0, α: $α, nmax: $nmax"

    # ∂c/∂t =  ∇ ⋅ (D ∇ c)
    # G(x,x',t) = det(D)^(-½) (4π t)^(-n/2)  exp( - (x -x')ᵀ D⁻¹ (x -x')ᵀ / (4t))

    # G(x,x',t) = det(D)^(-½) (4π t)^(-n/2)  exp( - (x -x')ᵀ D⁻¹ (x -x')ᵀ / (4t))

    @debug "sum(c0): $(sum(c0))"

    DIVAnd.diffusion!(ivol, nus, α, nmax, c0, c)

    @debug "sum(c): $(sum(c))"

    detD = prod(len_adjusted.^2)
    t = α * nmax
    c = c * sqrt((4π * t)^n / detD)

    @debug "range of c: $(extrema(c))"

    itp = LinearInterpolation(gridx,c,extrapolation_bc = NaN);
    ci = itp.(x...)

    weighti = 1 ./ ci
    clamp!(weighti, 0, 1)

    return weighti
end

Random.seed!(12343)

nobs = 1000
nobs = 10000
#nobs = 30000
#nobs = 100000
#nobs = 10000000

obslon = 4*rand(nobs).^2
obslat = 4*rand(nobs).^2

x = (obslon,obslat)

sel  = 1:10000
sel = 1:length(obslon)
len = (0.1,0.1)
x = (obslon,obslat)

weight = DIVAnd.weight_RtimesOne(x, len);
weighti = weight_RtimesOne_binning(x, len)

ratio = sqrt(mean((weighti-weight).^2)) / sqrt(mean(weighti.^2))
@debug "weight ratio: $ratio"


@test sqrt(mean((weighti-weight).^2)) < 0.3 * sqrt(mean(weighti.^2))
