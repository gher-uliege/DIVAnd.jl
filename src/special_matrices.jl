
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

