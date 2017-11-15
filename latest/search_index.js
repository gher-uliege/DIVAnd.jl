var documenterSearchIndex = {"docs": [

{
    "location": "index.html#",
    "page": "divand.jl documentation",
    "title": "divand.jl documentation",
    "category": "page",
    "text": "divand"
},

{
    "location": "index.html#divand.divandrun",
    "page": "divand.jl documentation",
    "title": "divand.divandrun",
    "category": "Function",
    "text": "divandrun(mask,pmn,xi,x,f,len,epsilon2; <keyword arguments>)\n\nPerform an n-dimensional variational analysis of the observations f located at the coordinates x. The array fi represent the interpolated field at the grid defined by the coordinates xi and the scales factors pmn.\n\nInput:\n\nmask: binary mask delimiting the domain. true is inside and false outside. For oceanographic application, this is the land-sea mask.\npmn: scale factor of the grid. pmn is a tuple with n elements. Every      element represents the scale factor of the corresponding dimension. Its      inverse is the local resolution of the grid in a particular dimension.\nxi: tuple with n elements. Every element represents a coordinate of the final grid on which the observations are interpolated\nx: tuple with n elements. Every element represents a coordinate of the observations\nf: value of the observations minus the background estimate (m-by-1 array).   (see note)\nlen: correlation length\nepsilon2: error variance of the observations (normalized by the error variance of the background field). epsilon2 can be a scalar (all observations have the same error variance and their errors are decorrelated), a vector (all observations can have a difference error variance and their errors are decorrelated) or a matrix (all observations can have a difference error variance and their errors can be correlated). If epsilon2 is a scalar, it is thus the inverse of the signal-to-noise ratio.\n\nOptional input arguments specified as keyword arguments\n\nvelocity: velocity of advection constraint. The default is      no-advection constraint\nalpha: alpha is vector of coefficients multiplying various terms in the      cost function. The first element multiplies the norm.      The other i-th element of alpha multiplies the (i+1)-th derivative.      Per default, the highest derivative is m = ceil(1+neff/2) where neff is the      effective dimension of the problem (the number of dimensions with a nonzero      correlation length).\n The values of alpha is the (m+1)th row of the Pascal triangle:\n    m=0         1\n    m=1       1   1\n    m=1     1   2   1     (n=1,2)\n    m=2   1   3   3   1   (n=3,4)\n    ...\nEOF, EOF: sub-space constraint. Orthogonal (EOF' WE^2 EOF = I) (units of      EOF: m^(-n/2))\nEOF_scaling, EOF_scaling: (dimensional)\nconstraints: a structure with user specified constrain\nmoddim: modulo for cyclic dimension (vector with n elements).    Zero is used for non-cyclic dimensions. Halo points should    not be included for cyclic dimensions. For example if the first dimension    is cyclic, then the grid point corresponding to mask(1,j) should be    between mask(end,1) (left neighbor) and mask(2,j) (right neighbor)\nfracindex: fractional indices (n-by-m array). If this array is specified,    then x and xi are not used.\ninversion: direct solver (:chol for Cholesky factorization) or a    interative solver (:pcg for preconditioned conjugate gradient) can be    used.\ncompPC: function that returns a preconditioner for the primal formulation    if inversion is set to 'pcg'. The function has the following arguments:\n     fun = compPC(iB,H,R)\nwhere iB is the inverse background error covariance, H the observation   operator and R the error covariance of the observation. The function compPC returns the   preconditioner fun(x,fx) computing fx = M  x (the inverse of M times x)   where M is a positive defined symmetric matrix [1].   Effectively, the system E⁻¹ A (E⁻¹)ᵀ (E x) = E⁻¹ b is solved for (E x) where E Eᵀ = M.   Ideally, M should this be similar to A, so that E⁻¹ A (E⁻¹)ᵀ is close to the identity matrix.\nfi0: starting field for iterative primal algorithm (same size as mask)\nf0: starting field for iterative dual algorithm (same size as the observations f)\noperatortype: Val{:sparse} for using sparse matrices (default) or Val{:MatFun} or using functions   to define the constrains.\nscale_len: true (default) if the correlation length-scale should be scaled such that the analysical   kernel reaches 0.6019072301972346 (besselk(1.,1.)) at the same distance. The kernel behaves thus similar to   the default kernel in two dimensions (alpha = [1,2,1]).\nalphabc : numerical value defining how the last grid points are stretched outward. 1, the default value mimics an infinite domain.\n\nTo have previous behaviour of finite domain use alphabc=0\n\nbtrunc : if provided defines where to truncate the calculation of the covariance matrix B. Only values up and including alpha[btrunc] will be calculated. IF the\n\n			iterative solution is calculated, the missing terms will be calculated on the fly during the conjugate gradient calulcations. Default value is none and full covariance calculation.\n\nOutput:\n\nfi: the analysed field\ns: structure with an array s.P representing the analysed error covariance\n\nNote:\n\nIf zero is not a valid first guess for your variable (as it is the case for   e.g. ocean temperature), you have to subtract the first guess from the   observations before calling divand and then add the first guess back in.\n\nExample:\n\nsee divand_simple_example.jl\n\nReferences\n\n[1]  https://en.wikipedia.org/w/index.php?title=Conjugate_gradient_method&oldid=761287292#The_preconditioned_conjugate_gradient_method\n\n\n\n"
},

{
    "location": "index.html#divand.jl-documentation-1",
    "page": "divand.jl documentation",
    "title": "divand.jl documentation",
    "category": "section",
    "text": "divandrun"
},

{
    "location": "index.html#Vocabulary-1",
    "page": "divand.jl documentation",
    "title": "Vocabulary",
    "category": "section",
    "text": "SDNCollection"
},

{
    "location": "index.html#Information-for-developpers-1",
    "page": "divand.jl documentation",
    "title": "Information for developpers",
    "category": "section",
    "text": ""
},

{
    "location": "index.html#Update-the-documentation-1",
    "page": "divand.jl documentation",
    "title": "Update the documentation",
    "category": "section",
    "text": "Installpip3 install --user mkdocs\npip3 install --user python-markdown-mathPkg.add(\"Documenter\")"
},

]}
