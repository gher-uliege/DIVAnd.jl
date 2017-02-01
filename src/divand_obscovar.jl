"""
R = divand_obscovar(epsilon2,m)

Create a matrix representing the observation error covariance R of size m x m. 

If epsilon2 is a scalar, then R = epsilon2 * I
If epsilon2 is a vector, then R = diag(epsilon2)
If epsilon2 is a matrix, then R = epsilon2
"""

function divand_obscovar(epsilon2,m)

if ndims(epsilon2) == 0
  R = Diagonal([epsilon2 for i=1:m]);
elseif ndims(epsilon2) == 1
  R = Diagonal(epsilon2)
else
  R = epsilon2
end

return R
end
