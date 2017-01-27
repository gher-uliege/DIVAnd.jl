"""
ind = statevector_sub2ind(sv,subscripts)

Compute from a tuple of subscripts the linear index in the packed state vector.
The first element of the subscript indicates the variable index and the remaining the spatial subscripts.
"""

function statevector_sub2ind(sv,subscripts)

# index of variable
ivar = subscripts[1]

# offset of variable
ioff = sv.ind[ivar]

# linear index in ivar-th array
index = sub2ind(sv.size[ivar],Tuple(subscripts[2:end])...)

return ioff + sv.unpacked2packed[ivar][index]

end
