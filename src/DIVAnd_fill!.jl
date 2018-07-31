# Adapted from
# from http://julialang.org/blog/2016/02/iteration
# to  fill values in a regular grid array.
# The major use to go from a coarse resolution to a fine resolution with grid ratio 2 i
# It should be quick and better than linear interpolation
# For the moment just weighted average using the nine points around ...
#

function DIVAnd_fill!(A::AbstractArray,B::AbstractArray,fillvalue)
    ntimes=1
    nd=ndims(A)

    # central weight
    cw=1

    RI =
        @static if VERSION >= v"0.7.0-beta.0"
            CartesianIndices(ntuple(i -> 1:size(A,i),nd))
        else
            CartesianRange(size(A))
        end

    I1, Iend = first(RI), last(RI)
    stencil = 3*one(CartesianIndex(I1))

    for nn=1:ntimes
        for indI in RI
            w = 0.0
            s = zero(eltype(A))

            B[indI] = A[indI]
            if !isequal(A[indI],fillvalue)

                RJ =
                    @static if VERSION >= v"0.7.0-beta.0"
                        # https://github.com/JuliaLang/julia/issues/15276#issuecomment-297596373
                        # let block work-around

                        let indI = indI, I1 = I1, I1 = I1, Iend = Iend, stencil = stencil
                            CartesianIndices(ntuple(
                                i-> max(I1[i], indI[i]-stencil[i]):min(Iend[i], indI[i]+stencil[i]),nd))
                        end
                    else
                        CartesianRange(max(I1, indI-stencil), min(Iend, indI+stencil))
                    end

                for indJ in RJ

                    if !isequal(A[indJ],fillvalue)
                        s += A[indJ]
                        if (indI==indJ)
                            w += cw
                        else
                            w += 1.
                        end
                    end
                    # end if not fill value
                end
                if w>0.0
                    B[indI] = s/w
                end
            end
        end

	    A .= B
    end

    return B
end
