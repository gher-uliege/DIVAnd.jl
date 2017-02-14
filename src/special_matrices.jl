
type CovarIS{T} <: AbstractMatrix{T}
    IS:: AbstractMatrix{T}
    factors
end

function CovarIS{T}(IS::AbstractMatrix{T})
    factors = nothing
    CovarIS(IS,factors)
end


Base.inv{T}(C::CovarIS{T}) = C.IS

Base.size{T}(C::CovarIS{T}) = size(C.IS)

function Base.:*{T}(C::CovarIS{T}, v::AbstractVector{Float64})
    if C.factors != nothing
        return C.factors \ v
    else
        return C.IS \ v
    end
end

Base.:*{T}(C::CovarIS{T}, v::SparseVector{Float64,Int64}) = C*full(v)


function A_mul_B{T}(C::CovarIS{T}, M::AbstractMatrix{Float64})
    if C.factors != nothing
        return C.factors \ M
    else
        return C.IS \ M
    end
end

# call to C * M
Base.:*{T}(C::CovarIS{T}, M::AbstractMatrix{Float64}) = A_mul_B(C,M)

# The following two definitions are necessary; otherwise the full C matrix will be formed when
# calculating C * M' or C * M.'

# call to C * M' (conjugate transpose: C Mᴴ)
Base.A_mul_Bc{T}(C::CovarIS{T}, M::AbstractMatrix{Float64}) = A_mul_B(C,M')
# call to C * M.' (transpose: C Mᵀ)
Base.A_mul_Bt{T}(C::CovarIS{T}, M::AbstractMatrix{Float64}) = A_mul_B(C,M.')


function Base.getindex{T}(C::CovarIS{T}, i::Int,j::Int)
    ei = zeros(eltype(C),size(C,1)); ei[i] = 1
    ej = zeros(eltype(C),size(C,1)); ej[j] = 1

    return (ej'*(C*ei))[1]
end


Base.:\{T}(C::CovarIS{T}, M::AbstractArray{Float64,2}) = C.IS * M

function factorize!{T}(C::CovarIS{T})
#    C.factors = cholfact(Symmetric(C.IS), Val{true})
    C.factors = cholfact(Symmetric(C.IS))
#    C.factors = cholfact(C.IS, Val{true})
end


function diagMtCM{T}(C::CovarIS{T}, M::AbstractMatrix{Float64})
    if C.factors != nothing
        return squeeze(sum((abs(C.factors[:PtL]\M)).^2,1),1)
    else
        return diag(M'*(C.IS \ M))
    end
end

function diagLtCM{T}(L::AbstractMatrix{Float64}, C::CovarIS{T}, M::AbstractMatrix{Float64})
    if C.factors != nothing
        return squeeze(sum((C.factors[:PtL]\M).*(C.factors[:PtL]\L),1),1)
    else
        return diag(L'*(C.IS \ M))
    end
end



# MatFun: a matrix defined by a function representing the matrix product

type MatFun{T}  <: AbstractMatrix{Float64}
    sz::Tuple{T,T}
    fun:: Function
    funt:: Function
end

Base.size{T}(MF::MatFun{T}) = MF.sz

for op in [:+, :-]; @eval begin
    function Base.$op{T}(MF1::MatFun{T}, MF2::MatFun{T})
        return MatFun(size(MF1),x -> $op(MF1.fun(x),MF2.fun(x)), x -> $op(MF2.funt(x),MF1.funt(x)))
    end

    Base.$op{T}(MF::MatFun{T}, S::AbstractSparseMatrix) = $op(MF,MatFun(S))
    Base.$op{T}(S::AbstractSparseMatrix, MF::MatFun{T}) = $op(MatFun(S),MF)
end
end

Base.:*{T,S}(MF::MatFun{T}, x::AbstractVector{S}) = MF.fun(x)
Base.:*{T,S}(MF::MatFun{T}, M::AbstractMatrix{S}) = cat(2,[MF.fun(M[:,i]) for i = 1:size(M,2)]...)

function A_mul_B{T}(MF1::MatFun{T}, MF2::MatFun{T})
    if size(MF1,2) != size(MF2,1)
        error("incompatible sizes")
    end
    return MatFun((size(MF1,1),size(MF2,2)),x -> MF1.fun(MF2.fun(x)), x -> MF2.funt(MF1.funt(x)))
end

Base.:*{T}(MF1::MatFun{T}, MF2::MatFun{T}) = A_mul_B(MF1,MF2)
Base.:*{T}(MF::MatFun{T}, S::AbstractSparseMatrix) = MF * MatFun(S)
Base.:*{T}(S::AbstractSparseMatrix,MF::MatFun{T}) = MatFun(S) * MF

for op in [:/, :*]; @eval begin
    Base.$op{T}(MF::MatFun{T}, a::Number) = MatFun(size(MF),x -> $op(MF.fun(x),a),x -> $op(MF.funt(x),a))
end
end

Base.:*{T}(a::Number, MF::MatFun{T}) = MatFun(size(MF),x -> a*MF.fun(x),x -> a*MF.funt(x))


function Base.:^{T}(MF::MatFun{T},n::Int)
    if n == 0
        return MatFun(size(MF),identity,identity)
    else
        return MF * (MF^(n-1))
    end
end

Base.:transpose{T}(MF:: MatFun{T}) = MatFun((MF.sz[2],MF.sz[1]),MF.funt,MF.fun)

Base.Ac_mul_B{T,S}(MF:: MatFun{T}, x::AbstractVector{S}) = MF.funt(x)
Base.Ac_mul_B{T}(MF1:: MatFun{T}, MF2:: MatFun{T}) = A_mul_B(MF1',MF2)


MatFun(S::AbstractSparseMatrix) = MatFun(size(S), x -> S*x, x -> S'*x)


# CovarHPHt representing H P Hᵀ

type CovarHPHt{T} <: AbstractMatrix{T}
    P:: AbstractMatrix{T}
    H:: AbstractMatrix{T}
end

Base.size{T}(C::CovarHPHt{T}) = (size(C.H,1),size(C.H,1))

#function CovarHPHt{T}(P::AbstractMatrix{T},H::AbstractMatrix{T})
#    CovarIS(IS,factors)
#end

function Base.:*{T}(C::CovarHPHt{T}, v::AbstractVector{Float64})
    return C.H * (C.P * (C.H' * v))
end


function A_mul_B{T}(C::CovarHPHt{T}, M::AbstractMatrix{Float64})
    return C.H * (C.P * (C.H' * M))
end

# call to C * M
Base.:*{T}(C::CovarHPHt{T}, M::AbstractMatrix{Float64}) = A_mul_B(C,M)

# The following two definitions are necessary; otherwise the full C matrix will be formed when
# calculating C * M' or C * M.'

# call to C * M' (conjugate transpose: C Mᴴ)
Base.A_mul_Bc{T}(C::CovarHPHt{T}, M::AbstractMatrix{Float64}) = A_mul_B(C,M')
# call to C * M.' (transpose: C Mᵀ)
Base.A_mul_Bt{T}(C::CovarHPHt{T}, M::AbstractMatrix{Float64}) = A_mul_B(C,M.')



function Base.getindex{T}(C::CovarHPHt{T}, i::Int,j::Int)
    ei = zeros(eltype(C),size(C,1)); ei[i] = 1
    ej = zeros(eltype(C),size(C,1)); ej[j] = 1

    return (ej'*(C*ei))[1]
end
