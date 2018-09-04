if VERSION >= v"0.7"
    using SparseArrays
    using LinearAlgebra
end



if VERSION >= v"0.7"
    Base.:*(A::SparseMatrixCSC{Float64,Int},B::Adjoint{Float64,SparseMatrixCSC{Float64,Int}}) = SparseArrays.spmatmul(A,copy(B),sortindices=:doubletranspose)
    Base.:*(A::Adjoint{Float64,SparseMatrixCSC{Float64,Int}},B::SparseMatrixCSC{Float64,Int}) = SparseArrays.spmatmul(copy(A),B,sortindices=:doubletranspose)
    Base.:*(A::SparseMatrixCSC{Float64,Int},B::SparseMatrixCSC{Float64,Int}) = SparseArrays.spmatmul(A,B,sortindices=:doubletranspose)
else
    Base.A_mul_Bc(A::SparseMatrixCSC{Float64,Int},B::SparseMatrixCSC{Float64,Int}) = SparseArrays.spmatmul(A,B',sortindices=:doubletranspose)
    Base.Ac_mul_B(A::SparseMatrixCSC{Float64,Int},B::SparseMatrixCSC{Float64,Int}) = SparseArrays.spmatmul(A',B,sortindices=:doubletranspose)
    Base.:*(A::SparseMatrixCSC{Float64,Int},B::SparseMatrixCSC{Float64,Int}) = SparseArrays.spmatmul(A,B,sortindices=:doubletranspose)
end
