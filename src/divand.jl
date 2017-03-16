module divand
using Interpolations
using NetCDF
using Roots
import SpecialFunctions
using Base.Test

type divand_constrain
    yo
    R
    H
end

type divand_struct
    n
    neff
    coeff
    sv
    D
    mask
    WE
    isinterior
    isinterior_stag
    isinterior_unpacked
    mapindex_packed
    mask_stag
    WEs
    WEss
    Dx
    alpha
    iB
    iB_
    moddim
    iscyclic
    applybc
    betap
    EOF_lambda
    primal
    factorize
    tol
    maxit
    minit
    inversion
    niter
    compPC
    progress
    preconditioner
    keepLanczosVectors
    yo
    R
    H
    P
    obsout
    obsconstrain

    function divand_struct(mask)
        n = ndims(mask)
        neff = 0
        coeff = 1.
        moddim = []
        iscyclic = []
        alpha = []
        yo = []
        R = []
        H = []

        sv = statevector_init((mask,))
        sz = size(mask)
        sempty = sparse(Array{Int64}([]),Array{Int64}([]),Array{Float64}([]),prod(sz),prod(sz))

        D = copy(sempty)
        WE = copy(sempty)
        WE = copy(sempty)
        iB = copy(sempty)
        iB_ = Array{SparseMatrixCSC{Float64,Int64}}(3);
        P = []

        isinterior = []
        isinterior_stag = [[] for i in 1:n]
        isinterior_unpacked = []
        mapindex_packed = []
        mask_stag = [Array{Bool,1}() for i in 1:n]
        WEs = [copy(sempty) for i in 1:n]
        WEss = [sparse(Array{Int64}([]),Array{Int64}([]),Array{Float64}([])) for i in 1:n]
        Dx = [copy(sempty) for i in 1:n]
        applybc = copy(sempty)

        betap = 0
        EOF_lambda = 0
        primal = true
        factorize = true
        tol = 1e-6
        maxit = 100
        minit = 10
        inversion = :chol
        niter = 0
        keepLanczosVectors = false

        obsout = Array{Bool,1}()
        obsconstrain = divand_constrain([],[],[])

        WEs = Array{Any,1}(n)
        WEss = Array{Any,1}(n)


        compPC(iB,R,H) = identity
        progress(iter,x,r,tol2,fun,b) = nothing
        preconditioner = identity

        new(n,
            neff,
            coeff,
            D,
            sv,
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
end

function ndgrid_fill(a, v, s, snext)
    for j = 1:length(a)
        a[j] = v[div(rem(j-1, snext), s)+1]
    end
end

# type stable
function ndgrid{T}(v1::AbstractVector{T},v2::AbstractVector{T})
    return ([x1 for x1 in v1, x2 in v2],
            [x2 for x1 in v1, x2 in v2])
end


function ndgrid{T}(vs::AbstractVector{T}...)
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
include("statevector.jl");
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
include("divand_Lpmnrange.jl");
include("divand_bc_stretch.jl");

include("divand_save.jl");

include("load_mask.jl");


export MatFun,divand_obscovar,divand_pc_sqrtiB,divand_pc_none,sparse_diag, statevector, pack, unpack, ind2sub, sub2ind, CovarHPHt, divand_rectdom, divand_squaredom, load_mask, oper_diag, oper_stagger, oper_diff, oper_pack, oper_trim, oper_shift, divand_save

export sparse_stagger, sparse_diff, localize_separable_grid, ndgrid, sparse_pack, sparse_interp, sparse_trim, sparse_shift, sparse_gradient, divand_laplacian,
statevector_init, statevector_pack, statevector_unpack, statevector_ind2sub, statevector_sub2ind, divandrun, divand_metric, distance, CovarIS, factorize!, divand_kernel, divand_cpme, divand_aexerr, divand_GCVKii, divand_diagHK, divand_GCVKiiobs, divand_diagHKobs, diagMtCM, diagLtCM, divand_residual, divand_residualobs,
divand_cvestimator, divand_erroratdatapoints, divand_cv, divand_qc, divand_adaptedeps2, divandgo, divandjog, divand_sampler, divand_filter3, divand_Lpmnrange, divand_cutter, divand_fittocpu, divand_bc_stretch

end
