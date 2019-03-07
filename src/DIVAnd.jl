VERSION < v"0.7.0-beta2.199" && __precompile__()


module DIVAnd

using Interpolations
using NCDatasets
using Base.Cartesian
using DataStructures
using SpecialFunctions
import HTTP
import Mustache
import ZipFile
using Missings
using StatsBase
#import JLD2, FileIO

if VERSION >= v"0.7.0-beta.0"
    using Printf
    using LinearAlgebra
    using SparseArrays
    using SuiteSparse
    using Distributed
    using Random
    using Statistics
    using SharedArrays
    using UUIDs
    using Dates
    using DelimitedFiles
    using Test

    # workaround for
    # https://github.com/JuliaLang/julia/issues/28011
    import Base: *
    Base.:*(A::SparseArrays.SparseMatrixCSC,B::BitArray) = A*Int8.(B)

    # surpress any logging info for testing
    macro nologs(expr)
        return :( @test_logs (:info,r".*") (:warn,r".*") match_mode=:any $expr )
    end
else
    using Base.Test
    const uuid1 = Base.Random.uuid1
    const mul! = A_mul_B!

    using Compat: @info, @warn, @debug
    using Compat: stdout, range, cat

    macro distributed(expr)
        return :( @parallel $(esc(expr)) )
    end
end

if (v"1.0" <= VERSION) && (VERSION < v"1.1")
    # performance regression in julia 1.0, but fixed in 1.1.0-DEV.308
    import Base:*
    (*)(A::SparseMatrixCSC, D::Diagonal) = A*sparse(D)
    (*)(D::Diagonal, A::SparseMatrixCSC) = sparse(D)*A
end


using Compat

const EarthRadius = 6372795.477598; # m

include("statevector.jl")
include("special_matrices.jl");

mutable struct DIVAnd_constrain{T <: AbstractFloat, TR <: AbstractMatrix{<: Number}, TH <: AbstractMatrix{<: Number}}
    yo::Vector{T}
    R::TR
    H::TH
end

# T is the type of floats and
# Ti: the type of integers
# N: the number of dimensions
mutable struct DIVAnd_struct{T,Ti,N,OT}
#mutable struct DIVAnd_struct{T,Ti,N,OT} where T <: AbstractFloat where Ti <: Int where N where OT
    n::Ti
    neff::Ti
    coeff::T
    sv::statevector{1,N}
    D::OT
    mask::BitArray{N}
    WE::OT
    isinterior::Vector{Bool}
    isinterior_stag::Vector{Vector{Bool}}
    isinterior_unpacked::Vector{Bool}
    mapindex_packed::Vector{Ti}
    mask_stag::Vector{BitArray{N}}
    WEs::Vector{OT}
    WEss::Vector{OT}
    Dx::NTuple{N,AbstractMatrix{T}}
    alpha::Vector{T}
    coeff_laplacian::Vector{T}
    coeff_derivative2::Vector{T}
    iB::OT
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
    R::Diagonal{T}
    #R::Union{Diagonal{T},SparseMatrixCSC{T,Int}}
    H::SparseMatrixCSC{T,Int}
    P::CovarIS{T,OT}
    obsout::BitArray{1}
    obsconstrain::DIVAnd_constrain{T}
end

function myempty(::Type{SparseMatrixCSC{T,Int}},sz) where T
    return sparse(Int[],Int[],Float64[],prod(sz),sz[1],sz[2])
end

function myempty(::Type{MatFun{Int}},sz)
    return MatFun(sz, x -> zero(x), x -> zero(x))
end


    function DIVAnd_struct(::Type{OT},mask::AbstractArray{Bool,n}) where {OT,n}
        neff = 0
        coeff = 1.
        moddim = Float64[]
        iscyclic = convert(Vector{Bool},falses(n))
        alpha = Float64[]
        coeff_laplacian = Float64[]
        coeff_derivative2 = Float64[]
        yo = Float64[]
        R = Diagonal(Float64[])
        H = sparse(Int[],Int[],Float64[],0,0)

        sv = statevector_init((mask,))
        sz = size(mask)
        sempty = myempty(OT,(prod(sz),prod(sz)))

        D = myempty(OT,(prod(sz),prod(sz)))
        WE = myempty(OT,(prod(sz),prod(sz)))
        iB = myempty(OT,(prod(sz),prod(sz)))
        iB_ = OT[]
        Ld = Float64[]
        P = CovarIS{Float64,OT}(sempty,nothing)

        isinterior = Bool[]
        isinterior_stag = [Bool[] for i in 1:n]
        isinterior_unpacked = Bool[]
        mapindex_packed = Int[]
        mask_stag = Vector{BitArray{n}}(undef,n)
        WEs = Array{OT,1}(undef,n)
        WEss = Array{OT,1}(undef,n)
        Dx = ([myempty(OT,(prod(sz),prod(sz))) for i in 1:n]...,)
        applybc = myempty(OT,(prod(sz),prod(sz)))

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

        obsout = BitArray{1}()
        obsconstrain = DIVAnd_constrain(Float64[],Matrix{Float64}(undef,0,0),Matrix{Float64}(undef,0,0))

        compPC(iB,R,H) = identity
        progress(iter,x,r,tol2,fun,b) = nothing
        preconditioner = identity

        return DIVAnd_struct{Float64,Int,n,OT}(n,
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
            coeff_laplacian,
            coeff_derivative2,
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

function DIVAnd_struct(mask::AbstractArray{Bool,n}) where n
    DIVAnd_struct(SparseMatrixCSC{Float64,Int},mask)
end

function DIVAnd_struct(::Type{Val{:sparse}},mask::AbstractArray{Bool,n}) where n
    DIVAnd_struct(SparseMatrixCSC{Float64,Int},mask)
end

function DIVAnd_struct(::Type{Val{:MatFun}},mask::AbstractArray{Bool,n}) where n
    DIVAnd_struct(MatFun{Int},mask)
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
    out = ntuple(i->Array{T}(undef,sz), n)
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
    Diagonal(vcat([diag(x) for x in X]...))
end

"""display size as a string """
formatsize(sz) = join(sz,"×")

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
function len_harmonize(len::T,mask::AbstractArray{Bool,N})::NTuple{N, Array{T,N}} where T <: Number where N
    return ((fill(len,size(mask)) for i=1:N)...,)
end

function len_harmonize(len::NTuple{N,T},mask::AbstractArray{Bool,N})::NTuple{N, Array{T,N}} where T <: Number where N
    return ((fill(len[i],size(mask)) for i=1:N)...,)
end

function len_harmonize(len::NTuple{N,AbstractArray{T,N}},mask::AbstractArray{Bool,N})::NTuple{N, Array{T,N}} where T <: Number where N
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

include("conjugategradient.jl");
include("DIVAnd_laplacian.jl");
include("DIVAnd_operators.jl");
include("DIVAnd_background_components.jl")
include("DIVAnd_background.jl");
include("DIVAnd_addc.jl");
include("DIVAnd_kernel.jl");
include("DIVAnd_obscovar.jl");
include("DIVAnd_obs.jl");
include("DIVAnd_pc_none.jl");
include("DIVAnd_pc_sqrtiB.jl");
include("DIVAnd_factorize.jl");
include("DIVAnd_solve.jl");
include("DIVAnd_metric.jl");
include("DIVAnd_constr_advec.jl");
include("DIVAnd_constr_fluxes.jl");
export DIVAnd_constr_fluxes

include("DIVAnd_constr_constcoast.jl");
export DIVAnd_constr_constcoast

include("DIVAndjog.jl");
include("DIVAndrun.jl");
include("DIVAndgo.jl");
include("DIVAnd_cpme.jl");
include("DIVAnd_aexerr.jl");
include("DIVAnd_GCVKii.jl");
include("DIVAnd_diagHK.jl");
include("DIVAnd_GCVKiiobs.jl");
include("DIVAnd_diagHKobs.jl");
include("DIVAnd_residual.jl");
include("DIVAnd_residualobs.jl");
include("DIVAnd_cvestimator.jl");
include("DIVAnd_erroratdatapoints.jl");
include("DIVAnd_cv.jl");
include("DIVAnd_qc.jl");
include("DIVAnd_sampler.jl");
include("DIVAnd_cutter.jl");
include("DIVAnd_fittocpu.jl");
include("DIVAnd_adaptedeps2.jl");
include("DIVAnd_filter3.jl");
include("DIVAnd_fill!.jl");
include("DIVAnd_Lpmnrange.jl");
include("DIVAnd_bc_stretch.jl");
include("DIVAnd_averaged_bg.jl")
include("jmBix.jl");

include("DIVAnd_iBpHtiRHx!.jl")

include("DIVAnd_save.jl");
export saveobs

include("varanalysis.jl");

include("load_mask.jl");
export extract_bath, load_bath, load_mask

include("load_obs.jl");
export loadbigfile, loadobs

include("domain.jl");
export domain, DIVAnd_rectdom, DIVAnd_squaredom

# high-level interface
include("diva.jl");
export diva3d

include("DIVAnd_weights.jl");

include("obsstat.jl");
export statpos

include("anamorphosis.jl");
export Anam

# ODV support
include("ODVspreadsheet.jl");
export ODVspreadsheet

# NetCDF support
include("NCSDN.jl");
export NCSDN

include("NCODV.jl");
export NCODV

# Vocabulary
include("Vocab.jl");
export Vocab
export urn_str

include("SDNMetadata.jl");
export SDNMetadata, SDNObsMetadata, divadoxml

include("select_time.jl");
export TimeSelectorYW, TimeSelectorYearListMonthList, TimeSelectorRunningAverage

include("fit.jl");
export fit_isotropic
export fitvertlen
export fithorzlen

include("utils.jl");

include("Quadtrees.jl");
export checkduplicates

include("DIVAnd_datainboundingbox.jl")

include("DIVAnd_cpme_go.jl")

include("scaleseparation.jl")


export DIVAnd_laplacian_prepare, DIVAnd_laplacian_apply, DIVAndrunfi

# statevector
export packens, unpackens

export MatFun,DIVAnd_obscovar,DIVAnd_pc_sqrtiB,DIVAnd_pc_none,sparse_diag, statevector, pack, unpack, ind2sub, sub2ind, CovarHPHt, oper_diag, oper_stagger, oper_diff, oper_pack, oper_trim, oper_shift, DIVAnd_save, varanalysis, dvmaskexpand, jmBix, DIVAnd_iBpHtiRHx!

export sparse_stagger, sparse_diff, localize_separable_grid, ndgrid, sparse_pack, sparse_interp, sparse_trim, sparse_shift, DIVAnd_gradient, sparse_gradient, DIVAnd_laplacian,
statevector_init, statevector_pack, statevector_unpack, statevector_ind2sub, statevector_sub2ind, DIVAndrun, DIVAnd_metric, distance, CovarIS, factorize!, DIVAnd_kernel, DIVAnd_cpme, DIVAnd_aexerr, DIVAnd_GCVKii, DIVAnd_diagHK, DIVAnd_GCVKiiobs, DIVAnd_diagHKobs, diagMtCM, diagLtCM, DIVAnd_residual, DIVAnd_residualobs,
DIVAnd_cvestimator, DIVAnd_erroratdatapoints, DIVAnd_cv, DIVAnd_qc, DIVAnd_adaptedeps2, DIVAndgo, DIVAndjog, DIVAnd_sampler, DIVAnd_filter3, DIVAnd_fill!, DIVAnd_Lpmnrange, DIVAnd_cutter, DIVAnd_fittocpu, DIVAnd_bc_stretch, DIVAnd_averaged_bg, DIVAnd_datainboundingbox, DIVAnd_cpme_go, scaleseparation

export checkobs


# old function name (to be depreciated)
const sparse_gradient = DIVAnd_gradient
export sparse_gradient
# @deprecate sparse_gradient DIVAnd_gradient

end
