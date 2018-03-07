#__precompile__()

module divand

using Interpolations
using NCDatasets
using DataArrays
using Base.Test
using Base.Cartesian
using DataStructures
import SpecialFunctions
import HTTP
using NLopt

include("statevector.jl")

type divand_constrain{T <: AbstractFloat, TR <: AbstractMatrix{<: Number}, TH <: AbstractMatrix{<: Number}}
    yo::Vector{T}
    R::TR
    H::TH
end

# T is the type of floats and
# Ti: the type of integers
# N: the number of dimensions
type divand_struct{T <: AbstractFloat,Ti <: Int,N}
    n::Ti
    neff::Ti
    coeff::T
    sv::statevector{1,N}
#    D::SparseMatrixCSC{T,Ti}
    D::AbstractMatrix{T}
    mask::BitArray{N}
    WE::AbstractMatrix{T}
    isinterior::Vector{Bool}
    isinterior_stag::Vector{Vector{Bool}}
    isinterior_unpacked::Vector{Bool}
    mapindex_packed::Vector{Ti}
    mask_stag::Vector{BitArray{N}}
    #WEs:: Vector{AbstractMatrix{T}} # does not work
    WEs::Vector{Any}
    WEss::Vector{Any}
    Dx::NTuple{N,AbstractMatrix{T}}
    alpha::Vector{T}
    iB::AbstractMatrix{T}
    #iB_::Vector{Any}
    iB_
    Ld::Array{T,1}
    moddim::Vector{T}
    iscyclic::Vector{Bool}
    applybc::AbstractMatrix{T}
    betap::T
    EOF_lambda::Vector{T}
    primal::Bool
    factorize::Bool
    tol::T
    maxit::Ti
    minit::Ti
    inversion::Symbol
    niter::Ti
    compPC
    progress
    preconditioner
    keepLanczosVectors::Bool
    yo::Vector{T}
    R::AbstractMatrix{T}
    H::AbstractMatrix{T}
    P::AbstractMatrix{T}
    obsout::BitArray{1}
    obsconstrain::divand_constrain{T}
end


    function divand_struct(mask)
        n = ndims(mask)
        neff = 0
        coeff = 1.
        moddim = Float64[]
        iscyclic = convert(Vector{Bool},falses(n))
        alpha = Float64[]
        yo = Float64[]
        R = Matrix{Float64}(0,0)
        H = Matrix{Float64}(0,0)

        sv = statevector_init((mask,))
        sz = size(mask)
        sempty = sparse(Array{Int}([]),Array{Int}([]),Array{Float64}([]),prod(sz),prod(sz))

        D = copy(sempty)
        WE = copy(sempty)
        iB = copy(sempty)
        iB_ = Vector{SparseMatrixCSC{Float64,Int}}()
        Ld = Float64[]
        P = Matrix{Float64}(0,0)

        isinterior = Bool[]
        isinterior_stag = [Bool[] for i in 1:n]
        isinterior_unpacked = Bool[]
        mapindex_packed = Int[]
        #mask_stag = [Array{Bool,1}() for i in 1:n]
        mask_stag = [BitArray{n}(zeros(Int,n)...) for i in 1:n]
        WEs = Vector{SparseMatrixCSC{Float64,Int}}()
        WEss = [sparse(Array{Int}([]),Array{Int}([]),Array{Float64}([])) for i in 1:n]
        Dx = ([sparse(Array{Int}([]),Array{Int}([]),Array{Float64}([])) for i in 1:n]...)
        applybc = copy(sempty)

        betap = 0.
        EOF_lambda = Float64[]
        primal = true
        factorize = true
        tol = 1e-6
        maxit = 100
        minit = 10
        inversion = :chol
        niter = 0
        keepLanczosVectors = false

        #obsout = Array{Bool,1}()
        obsout = BitArray{1}()
        obsconstrain = divand_constrain(Float64[],Matrix{Float64}(0,0),Matrix{Float64}(0,0))

        WEs = Array{Any,1}(n)
        WEss = Array{Any,1}(n)


        compPC(iB,R,H) = identity
        progress(iter,x,r,tol2,fun,b) = nothing
        preconditioner = identity

        return divand_struct(n,
            neff,
            coeff,
            sv,
            D,
            mask,
            WE,
            isinterior,
            isinterior_stag,
            isinterior_unpacked,
            mapindex_packed,
            mask_stag,
            WEs,
            WEss,
            Dx,
            alpha,
            iB,
            iB_,
            Ld,
            moddim,
            iscyclic,
            applybc,
            betap,
            EOF_lambda,
            primal,
            factorize,
            tol,
            maxit,
            minit,
            inversion,
            niter,
            compPC,
            progress,
            preconditioner,
            keepLanczosVectors,
            yo,
            R,
            H,
            P,
            obsout,
            obsconstrain
            )
    end

# ndgrid* functions are from
# https://github.com/JuliaLang/julia/blob/master/examples/ndgrid.jl
# Licence MIT

function ndgrid_fill(a, v, s, snext)
    for j = 1:length(a)
        a[j] = v[div(rem(j-1, snext), s)+1]
    end
end

# type stable
function ndgrid(v1::AbstractVector{T},v2::AbstractVector{T}) where T
    return ([x1 for x1 in v1, x2 in v2],
            [x2 for x1 in v1, x2 in v2])
end


function ndgrid(vs::AbstractVector{T}...) where T
    n = length(vs)
    sz = map(length, vs)
    out = ntuple(i->Array{T}(sz), n)
    s = 1
    for i=1:n
        a = out[i]::Array
        v = vs[i]
        snext = s*size(a,i)
        ndgrid_fill(a, v, s, snext)
        s = snext
    end
    out
end


# https://stackoverflow.com/questions/31235469/array-type-promotion-in-julia
function promote_array(arrays...)
    eltype = Base.promote_eltype(arrays...)
    tuple([convert(Array{eltype}, array) for array in arrays]...)
end

ndgrid(vs...) = ndgrid(promote_array(vs...)...)

"""concatenate diagonal matrices"""
function blkdiag(X::Diagonal...)
    Diagonal(cat(1,[diag(x) for x in X]...))
end

"""display size as a string """
formatsize(sz) = join(sz,"×")

"""
mask,xyi,pmn = divand_squaredom(n,coord)

Create a "square" domain in `n` dimensions with the coordinates `coord`
assuming a Catersian metric. This functions returns
the mask `mask`, the coordinates `(xi,yi,...)` and the metric `(pm,pn...)`.

# Example

mask,(pm,pn),(xi,yi) = divand_squaredom(2,linspace(0,1,50))
"""
function divand_squaredom(n,coord)
    coords = ([coord for i = 1:n]...)
    return divand_rectdom(coords...)
end


"""
mask,xyi,pmn = divand_squaredom(n,coord)

Create a "square" domain in `n` dimensions with the coordinates `coord`
assuming a Catersian metric. This functions returns
the mask `mask`, the coordinates `(xi,yi,...)` and the metric `(pm,pn...)`.

# Example

mask,(pm,pn),(xi,yi) = divand_rectdom(linspace(0,1,50),linspace(0,1,50))
"""
function divand_rectdom(coords...)
    # grid of background field
    xyi = ndgrid(coords...)

    # mask (all points are valid)
    mask = trues(xyi[1])

    # metric (inverse of the resolution)
    pmn = ([ones(size(mask)) / (coords[i][2]-coords[i][1]) for i = 1:length(coords)]...)

    return mask,pmn,xyi
end


function dvmaskexpand(x)
    z=deepcopy(x)
    sz=size(z)[1]
    for i=1:sz
        if !z[i]
            ip=min(i+1,sz)
            im=max(i-1,1)
            z[i]=(x[im] | x[ip])
        end
    end
    return z
end


"""
Len = len_harmonise(len,mask)
Produce a tuple of arrays of the correlation length `len` which can be either a scalar (homogeneous and isotropic case),
a tuple of scalar (homogeneous case) or already a tuple of arrays (general case). The the later case the size of the arrays are veryfied.
"""
function len_harmonize{T <: Number,N}(len::T,mask::AbstractArray{Bool,N})::NTuple{N, Array{T,N}}
    return ((fill(len,size(mask)) for i=1:N)...)
end

function len_harmonize{T <: Number,N}(len::NTuple{N,T},mask::AbstractArray{Bool,N})::NTuple{N, Array{T,N}}
    return ((fill(len[i],size(mask)) for i=1:N)...)
end

function len_harmonize{T <: Number,N}(len::NTuple{N,AbstractArray{T,N}},mask::AbstractArray{Bool,N})::NTuple{N, Array{T,N}}
    for i=1:N
        if size(mask) != size(len[i])
            error("mask ($(formatsize(size(mask)))) and correlation length ($(formatsize(size(len[i])))) have incompatible size")
        end
    end

    return len
end

function len_harmonize(len,mask::AbstractArray{Bool,N}) where N
    # promote all lens to a common type
    return len_harmonize(promote_array(len...),mask)
end

@inline function alpha_default(neff::Int)
    # kernel should has be continuous derivative

    # highest derivative in cost function
    m = ceil(Int,1+neff/2)

    # alpha is the (m+1)th row of the Pascal triangle:
    # m=0         1
    # m=1       1   1
    # m=1     1   2   1
    # m=2   1   3   3   1
    # ...

    return Int[binomial(m,k) for k = 0:m]
end


"""
    neff, alpha = alpha_default(Labs,alpha)

Return a default value of alpha.
"""

@inline function alpha_default(Labs,alpha::Vector{T}) where T
    # must handle the case when Labs is zero in some dimension
    # thus reducing the effective dimension
    neff = sum([mean(L) > 0 for L in Labs])::Int

    if isempty(alpha)
        return neff, Vector{T}(alpha_default(neff))
    else
        return neff, alpha
    end

end


include("sparse_operator.jl");
include("function_operator.jl");

include("sparse_stagger.jl");
include("sparse_diff.jl");
include("sparse_interp.jl");
include("sparse_trim.jl");
include("sparse_shift.jl");
include("sparse_gradient.jl");


include("localize_separable_grid.jl");

include("special_matrices.jl");
include("conjugategradient.jl");
include("divand_laplacian.jl");
include("divand_operators.jl");
include("divand_background_components.jl")
include("divand_background.jl");
include("divand_addc.jl");
include("divand_kernel.jl");
include("divand_obscovar.jl");
include("divand_obs.jl");
include("divand_pc_none.jl");
include("divand_pc_sqrtiB.jl");
include("divand_factorize.jl");
include("divand_solve.jl");
include("divand_metric.jl");
include("divand_constr_advec.jl");
include("divand_constr_fluxes.jl");
include("divandjog.jl");
include("divandrun.jl");
include("divandgo.jl");
include("divand_cpme.jl");
include("divand_aexerr.jl");
include("divand_GCVKii.jl");
include("divand_diagHK.jl");
include("divand_GCVKiiobs.jl");
include("divand_diagHKobs.jl");
include("divand_residual.jl");
include("divand_residualobs.jl");
include("divand_cvestimator.jl");
include("divand_erroratdatapoints.jl");
include("divand_cv.jl");
include("divand_qc.jl");
include("divand_sampler.jl");
include("divand_cutter.jl");
include("divand_fittocpu.jl");
include("divand_adaptedeps2.jl");
include("divand_filter3.jl");
include("divand_fill!.jl");
include("divand_Lpmnrange.jl");
include("divand_bc_stretch.jl");
include("divand_averaged_bg.jl")
include("jmBix.jl");

include("divand_iBpHtiRHx!.jl")

include("divand_save.jl");

include("varanalysis.jl");

include("load_mask.jl");
include("load_obs.jl");
export loadbigfile

include("domain.jl");
export domain

# high-level interface
include("diva.jl");
export diva, diva3d

include("divand_weights.jl");

include("obsstat.jl");
export statpos

include("anamorphosis.jl");
export Anam

# ODV support
include("ODVspreadsheet.jl");
export ODVspreadsheet

# Vocabulary
include("Vocab.jl");
export Vocab
export urn_str

include("SDNMetadata.jl");
export SDNMetadata

include("select_time.jl");

include("fit.jl");



export divand_laplacian_prepare, divand_laplacian_apply, divandrunfi

# statevector
export packens, unpackens

export MatFun,divand_obscovar,divand_pc_sqrtiB,divand_pc_none,sparse_diag, statevector, pack, unpack, ind2sub, sub2ind, CovarHPHt, divand_rectdom, divand_squaredom, load_mask, oper_diag, oper_stagger, oper_diff, oper_pack, oper_trim, oper_shift, divand_save, varanalysis, dvmaskexpand, jmBix, divand_iBpHtiRHx!

export sparse_stagger, sparse_diff, localize_separable_grid, ndgrid, sparse_pack, sparse_interp, sparse_trim, sparse_shift, sparse_gradient, divand_laplacian,
statevector_init, statevector_pack, statevector_unpack, statevector_ind2sub, statevector_sub2ind, divandrun, divand_metric, distance, CovarIS, factorize!, divand_kernel, divand_cpme, divand_aexerr, divand_GCVKii, divand_diagHK, divand_GCVKiiobs, divand_diagHKobs, diagMtCM, diagLtCM, divand_residual, divand_residualobs,
divand_cvestimator, divand_erroratdatapoints, divand_cv, divand_qc, divand_adaptedeps2, divandgo, divandjog, divand_sampler, divand_filter3, divand_fill!, divand_Lpmnrange, divand_cutter, divand_fittocpu, divand_bc_stretch, divand_averaged_bg

end
