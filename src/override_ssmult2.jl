# adapted from Julia


# See also
# https://github.com/PetterS/SuiteSparse/blob/27e5a8516464a6ac40bd3fa0e5b46e51b11f4765/CHOLMOD/MatrixOps/cholmod_ssmult.c#L239

"""compute non-zero values of A*B"""
function myspmatmul_nnz(A::SparseMatrixCSC{Tv,Ti}, B::SparseMatrixCSC{Tv,Ti}) where {Tv,Ti}

    mA, nA = size(A)
    mB, nB = size(B)
    nA==mB || throw(DimensionMismatch())

    colptrA = A.colptr; rowvalA = A.rowval; nzvalA = A.nzval
    colptrB = B.colptr; rowvalB = B.rowval; nzvalB = B.nzval

    ip = 1::Ti

    @inbounds begin
        xb = zeros(Ti, mA)
        for i in 1:nB

            for jp in colptrB[i]:(colptrB[i+1] - 1)

                j = rowvalB[jp]
                for kp in colptrA[j]:(colptrA[j+1] - 1)
                    k = rowvalA[kp]
                    if xb[k] != i
                        ip += 1
                        xb[k] = i
                    end
                end
            end
        end
    end

    return ip-1
end



function myspmatmul_unsorted(A::SparseMatrixCSC{Tv,Ti}, B::SparseMatrixCSC{Tv,Ti}, nnzC::Ti) where {Tv,Ti}

    mA, nA = size(A)
    mB, nB = size(B)
    nA==mB || throw(DimensionMismatch())

    colptrA = A.colptr; rowvalA = A.rowval; nzvalA = A.nzval
    colptrB = B.colptr; rowvalB = B.rowval; nzvalB = B.nzval

    colptrC = Array{Ti}(nB+1)
    rowvalC = Array{Ti}(nnzC)
    nzvalC = Array{Tv}(nnzC)

    @inbounds begin
        ip = 1
        xb = zeros(Ti, mA)
        x  = zeros(Tv, mA)
        for i in 1:nB
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

    Cunsorted = SparseMatrixCSC(mA, nB, colptrC, rowvalC, nzvalC)
    return Cunsorted
end



function Base.:*(A::SparseMatrixCSC{Float64,Int},B::SparseMatrixCSC{Float64,Int})

    nnzC = myspmatmul_nnz(A,B)
    #@show nnz(A),nnz(B),nnzC
    if nnz(A) + nnz(B) < nnzC
        # compute A*B as (B'*A')'
        # swap and transpose        
        #@show "swap"
        return myspmatmul_unsorted(B',A',nnzC)'        
    else
        Cunsorted = myspmatmul_unsorted(A,B,nnzC)
        C = SparseArrays.sortSparseMatrixCSC!(Cunsorted, sortindices=:doubletranspose)
        return C
    end
end


function Base.:A_mul_Bc(A::SparseMatrixCSC{Float64,Int},B::SparseMatrixCSC{Float64,Int})

    BT = B'
    nnzC = myspmatmul_nnz(A,BT)
    #@show nnz(A),nnz(B),nnzC
    if nnz(A) + nnz(B) < nnzC
        # compute A*B' as (B*A')'
        # swap and transpose        
        #@show "swap"
        return myspmatmul_unsorted(B,A',nnzC)'        
    else
        Cunsorted = myspmatmul_unsorted(A,BT,nnzC)
        C = SparseArrays.sortSparseMatrixCSC!(Cunsorted, sortindices=:doubletranspose)
        return C
    end
end


function Base.:Ac_mul_B(A::SparseMatrixCSC{Float64,Int},B::SparseMatrixCSC{Float64,Int})

    AT = A'
    nnzC = myspmatmul_nnz(AT,B)
    #@show nnz(A),nnz(B),nnzC
    if nnz(A) + nnz(B) < nnzC
        # compute A'*B as (B'*A)'
        # swap and transpose        
        #@show "swap"
        return myspmatmul_unsorted(B',A,nnzC)'        
    else
        Cunsorted = myspmatmul_unsorted(AT,B,nnzC)
        C = SparseArrays.sortSparseMatrixCSC!(Cunsorted, sortindices=:doubletranspose)
        return C
    end
end


