module divand
using Interpolations
using Base.Test

function ndgrid_fill(a, v, s, snext)
    for j = 1:length(a)
        a[j] = v[div(rem(j-1, snext), s)+1]
    end
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


sparse_diag(d) = spdiagm(d)


function sparse_pack(mask)

j = find(mask)
m = length(j)
i = collect(1:m)
s = ones(m)
n = length(mask)
H = sparse(i,j,s,m,n)

end

include("sparse_stagger.jl");
include("sparse_diff.jl");
include("sparse_interp.jl");
include("sparse_trim.jl");
include("sparse_shift.jl");
include("sparse_gradient.jl");
include("divand_laplacian.jl");
include("localize_separable_grid.jl");


include("statevector_init.jl");
include("statevector_pack.jl");
include("statevector_unpack.jl");



function test()
    include("test_sparse_diff.jl");
    include("test_localize_separable_grid.jl");
    include("test_statevector.jl");


end

export test, sparse_stagger, sparse_diff, localize_separable_grid, ndgrid, sparse_pack, sparse_interp, sparse_trim, sparse_shift, sparse_gradient, divand_laplacian,
   statevector_init, statevector_pack, statevector_unpack

end
