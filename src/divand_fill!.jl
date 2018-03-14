# Adapted from
# from http://julialang.org/blog/2016/02/iteration
# to  fill values in a regular grid array. 
# The major use to go from a coarse resolution to a fine resolution with grid ratio 2 i
# It should be quick and better than linear interpolation 
# For the moment just weighted average using the nine points around ...
#

function divand_fill!(A::AbstractArray,B::AbstractArray,fillvalue)

    #
    function dvisvalue(x)
        if isnan(fillvalue)
            return !isnan(x);
        else
            return !(x==fillvalue);
        end
    end
    ntimes=1
    nd=ndims(A)
    # central weight
    cw=3^nd-1
    cw=1
    

    R = CartesianRange(size(A))
    I1, Iend = first(R), last(R)
    for nn=1:ntimes

        for I in R
            w, s = 0.0, zero(eltype(A))
            
            B[I] = A[I]
            if !dvisvalue(A[I])
                for J in CartesianRange(max(I1, I-3*I1), min(Iend, I+3*I1))
                    #@show J
                    if dvisvalue(A[J])
                        s += A[J]
                        if (I==J)
                            w += cw
                        else
                            w += 1.
                        end
                    end
                    # end if not fill value
                end
                if w>0.0
                    B[I] = s/w
                end
            end
        end
        
		
	A[:]=B[:]	
		
    end


    return B
end
