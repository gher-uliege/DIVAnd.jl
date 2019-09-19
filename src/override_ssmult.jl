using SparseArrays
using LinearAlgebra


Base.:*(A::SparseMatrixCSC{Float64,Int}, B::Adjoint{Float64,SparseMatrixCSC{Float64,Int}}) =
    SparseArrays.spmatmul(A, copy(B), sortindices = :doubletranspose)
Base.:*(A::Adjoint{Float64,SparseMatrixCSC{Float64,Int}}, B::SparseMatrixCSC{Float64,Int}) =
    SparseArrays.spmatmul(copy(A), B, sortindices = :doubletranspose)
Base.:*(A::SparseMatrixCSC{Float64,Int}, B::SparseMatrixCSC{Float64,Int}) =
    SparseArrays.spmatmul(A, B, sortindices = :doubletranspose)
