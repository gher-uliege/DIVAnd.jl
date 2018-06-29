# adapter from Julia

function myspmatmul(A::SparseMatrixCSC{Tv,Ti}, B::SparseMatrixCSC{Tv,Ti};
                    sortindices::Symbol = :sortcols) where {Tv,Ti}

    mA, nA = size(A)
    mB, nB = size(B)
    nA==mB || throw(DimensionMismatch())

    colptrA = A.colptr; rowvalA = A.rowval; nzvalA = A.nzval
    colptrB = B.colptr; rowvalB = B.rowval; nzvalB = B.nzval
    # TODO: Need better estimation of result space
    nnzC = min(mA*nB, length(nzvalA) + length(nzvalB))
    colptrC = Array{Ti}(nB+1)
    rowvalC = Array{Ti}(nnzC)
    nzvalC = Array{Tv}(nnzC)

    @inbounds begin
        ip = 1
        xb = zeros(Ti, mA)
        x  = zeros(Tv, mA)
        for i in 1:nB
            if ip + mA - 1 > nnzC
                resize!(rowvalC, nnzC + max(nnzC,mA))
                resize!(nzvalC, nnzC + max(nnzC,mA))
                nnzC = length(nzvalC)
            end
            colptrC[i] = ip
            for jp in colptrB[i]:(colptrB[i+1] - 1)
                nzB = nzvalB[jp]
                j = rowvalB[jp]
                for kp in colptrA[j]:(colptrA[j+1] - 1)
                    nzC = nzvalA[kp] * nzB
                    k = rowvalA[kp]
                    if xb[k] != i
                        rowvalC[ip] = k
                        ip += 1
                        xb[k] = i
                        x[k] = nzC
                    else
                        x[k] += nzC
                    end
                end
            end
            for vp in colptrC[i]:(ip - 1)
                nzvalC[vp] = x[rowvalC[vp]]
            end
        end
        colptrC[nB+1] = ip
    end

    deleteat!(rowvalC, colptrC[end]:length(rowvalC))
    deleteat!(nzvalC, colptrC[end]:length(nzvalC))

    # The Gustavson algorithm does not guarantee the product to have sorted row indices.
    Cunsorted = SparseMatrixCSC(mA, nB, colptrC, rowvalC, nzvalC)
    #C = SparseArrays.sortSparseMatrixCSC!(Cunsorted, sortindices=sortindices)
    C = SparseArrays.sortSparseMatrixCSC!(Cunsorted, sortindices=:doubletranspose)
    return C
end


Base.A_mul_Bc(A::SparseMatrixCSC{Float64,Int},B::SparseMatrixCSC{Float64,Int}) = myspmatmul(A,B')
Base.Ac_mul_B(A::SparseMatrixCSC{Float64,Int},B::SparseMatrixCSC{Float64,Int}) = myspmatmul(A',B)
Base.:*(A::SparseMatrixCSC{Float64,Int},B::SparseMatrixCSC{Float64,Int}) = myspmatmul(A,B)
