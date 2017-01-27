"""
subscripts = statevector_ind2ind(sv,index)

Compute from linear index in the packed state vector a tuple of subscripts.
The first element of the subscript indicates the variable index and the remaining the spatial subscripts.
"""

function statevector_ind2sub(sv,index)

# variable index
ivar = sum(sv.ind .< index)

# substract offset
vind = index - sv.ind[ivar]

# spatial subscript
subscript = ind2sub(sv.size[ivar],sv.packed2unpacked[ivar][vind])

return (ivar,subscript...)

end
