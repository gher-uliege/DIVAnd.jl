
sparse_diag(d)::SparseMatrixCSC{Float64,Int64} = spdiagm(d)

function sparse_pack(mask)

    j = find(mask)
    m = length(j)
    i = collect(1:m)
    s = ones(m)
    n = length(mask)
    H = sparse(i,j,s,m,n)
end
