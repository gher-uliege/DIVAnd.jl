
#=

The gradient in drection x at the point > should be close to zero if one 
of the grid cell denoted by * is a land point.

  +-----------------+-----------------+-----------------+-----------------+
  |                 |                 |                 |                 |
  |                 |                 |                 |                 |
  |                 |                 |                 |                 |
  |                 |        *        |        *        |                 |
  |                 |     (i,j+1)     |    (i+1,j+1)    |                 |
  |                 |                 |                 |                 |
  |                 |                 |                 |                 |
  +-----------------+--------^--------+-----------------+-----------------+
  |                 |      (i,j)      |                 |                 |
  |                 |                 |                 |                 |
  |                 |                 |                 |                 |
  |                 |        .        >                 |                 |
  |                 |      (i,j)    (i,j)               |                 |
  |                 |                 |                 |                 |
  |                 |                 |                 |                 |
  +-----------------+-----------------+-----------------+-----------------+
  |                 |                 |                 |                 |
  |                 |                 |                 |                 |
  |                 |                 |                 |                 |
  |                 |        *        |        *        |                 |
  |                 |     (i,j-1)     |    (i+1,j-1)    |                 |
  |                 |                 |                 |                 |
  |                 |                 |                 |                 |
  +-----------------+-----------------+-----------------+-----------------+

=#

function boundaryu(mask)

    bu = falses(size(mask,1)-1,size(mask,2))
    
    for j = 1:size(bu,2)
        for i = 1:size(bu,1)
            # check if i+½,j is a see point
            
            if mask[i,j] && mask[i+1,j]
            
                if j > 1
                    if !mask[i,j-1] || !mask[i+1,j-1]
                        bu[i,j] = true
                    end
                end
            
                if j < size(bu,2)
                    if !mask[i,j+1] || !mask[i+1,j+1]
                        bu[i,j] = true
                    end            
                end
            end
        end
        
    end
    return bu
end

boundaryv(mask) = boundaryu(mask')'


function gradcoast(mask,fip)
    fi = zeros(size(mask))
    fi[mask] = fip
    bu = boundaryu(mask)
    bv = boundaryv(mask)
    Dxfi = fi[2:end,:] - fi[1:end-1,:]
    Dyfi = fi[:,2:end] - fi[:,1:end-1]

    return [Dxfi[bu]; Dyfi[bv]]
end


"""
    H = sparse_gradcoast(mask)

Return a sparse matrix `H` computing the
gradient of the values along the coastline defined by `mask`.
Land points correspond to the value `false` and sea points to the value `true`
in `mask`.
"""
function sparse_gradcoast(mask)
    sz = size(mask)
    Hunpack = sparse_pack(mask)'
    bu = boundaryu(mask)
    bv = boundaryv(mask)

    diffx = sparse_diff(sz,1)
    diffy = sparse_diff(sz,2)
    
    H = [sparse_pack(bu) * diffx * Hunpack;
         sparse_pack(bv) * diffy * Hunpack]
    
    return H
end

"""
    c = DIVAnd_constr_constcoast(mask,eps2)

Constrain imposing that the gradients along the coastline defined
by `mask` are close to zero constrolled by the parameter `eps2` which
represents the scalled error variance on the gradients.

This constrain is useful to indirectly impose that a stream function
does not have a current component perpendicular to the coast line.
"""
function DIVAnd_constr_constcoast(mask,eps2)
    H = sparse_gradcoast(mask)
    m = size(H,1)
    yo = zeros(m)
    R = Diagonal(fill(eps2,(m,)))
    
    return DIVAnd_constrain(yo,R,H)
end


