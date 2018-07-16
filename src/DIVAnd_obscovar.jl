"""
R = DIVAnd_obscovar(epsilon2,m)

Create a matrix representing the observation error covariance R of size m x m.

If epsilon2 is a scalar, then R = epsilon2 * I
If epsilon2 is a vector, then R = diag(epsilon2)
If epsilon2 is a matrix, then R = epsilon2
"""
DIVAnd_obscovar(epsilon2::Number,m) = Diagonal(fill(epsilon2,m))
DIVAnd_obscovar(epsilon2::Vector,m) = Diagonal(epsilon2)
DIVAnd_obscovar(epsilon2::AbstractMatrix,m) = epsilon2

