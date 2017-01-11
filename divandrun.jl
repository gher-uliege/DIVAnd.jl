# Compute a variational analysis of arbitrarily located observations.
#
# [fi,err,s] = divand(mask,pmn,xi,x,f,len,lambda,...);
# [fi,err] = divand(mask,pmn,xi,x,f,len,lambda,...);
# [fi,s] = divand(mask,pmn,xi,x,f,len,lambda,...);
#
# Perform an n-dimensional variational analysis of the observations f located at
# the coordinates x. The array fi represent the interpolated field at the grid
# defined by the coordinates xi and the scales factors pmn.
#
# Input:
#   mask: binary mask delimiting the domain. 1 is inside and 0 outside.
#         For oceanographic application, this is the land-sea mask.
#
#   pmn: scale factor of the grid. pmn is a cell array with n elements. Every 
#        element represents the scale factor of the corresponding dimension. Its
#        inverse is the local resolution of the grid in a particular dimension.
#
#   xi: cell array with n elements. Every element represents a coordinate
#   of the final grid on which the observations are interpolated
#   x: cell array with n elements. Every element represents a coordinate of
#   the observations
#   f: value of the observations *minus* the background estimate (m-by-1 array).
#     (see note)
#   len: correlation length
#   lambda: signal-to-noise ratio of observations (if lambda is a scalar).
#     The larger this value is, the closer is the field fi to the
#     observation.
#     if lambda is a scalar:
#        R = 1/lambda I, where R is the observation error covariance
#        matrix),
#     if lambda is a vector
#     a vector (R = diag(lambda)) or a matrix or a CovarParam object (R = 
#     lambda).
#
#
# Optional input arguments specified as pairs of keyword and values:
#  'velocity', vel: velocity of advection constraint. The default is 
#        no-advection constraint
#
#  'alpha': alpha is vector of coefficients multiplying various terms in the 
#        cost function. The first element multiplies the norm.
#        The other i-th element of alpha multiplies the (i+1)-th derivative. 
#        Per default, the highest derivative is m = ceil(1+n/2) where n is the 
#        dimension of the problem.
#
#        The values of alpha is the (m+1)th row of the Pascal triangle:
#           m=0         1
#           m=1       1   1
#           m=1     1   2   1     (n=1,2)
#           m=2   1   3   3   1   (n=3,4)
#           ...
#
#  'diagnostics': 0 or 1 turns diagnostic and debugging information on (1) or 
#        off (0, default). If on, they will be returned as the last output 
#        argument
#
#  'EOF', EOF: sub-space constraint. Orthogonal (EOF' WE^2 EOF = I) (units of 
#        EOF: m^(-n/2))
#
#  'EOF_scaling', EOF_scaling: (dimensional)
#
#  'constraint': a structure with user specified constrain
#
#  'moddim': modulo for cyclic dimension (vector with n elements). 
#      Zero is used for non-cyclic dimensions. Halo points should 
#      not be included for cyclic dimensions. For example if the first dimension
#      is cyclic, then the grid point corresponding to mask(1,j) should be 
#      between mask(end,1) (left neighbor) and mask(2,j) (right neighbor)
#
#  'fracdim': fractional indices (n-by-m array). If this array is specified, 
#      then x and xi are not used.
#
#  'inversion': direct solver ('chol' for Cholesky factorization) or a 
#      interative solver ('pcg' for preconditioned conjugate gradient) can be 
#      used.
#
#  'compPC': function that returns a preconditioner for the primal formulation 
#      if inversion is set to 'pcg'. The function has the following arguments:
#
#            [M1,M2] = compPC(iB,H,R)
#
#     where iB is the inverse background error covariance, H the observation 
#     operator and R the error covariance of the observation. The used 
#     preconditioner M will be M = M1 * M2 (or just M = M1 if M2 is empty).
#     Per default a modified incomplete Cholesky factorization will be used a 
#     preconditioner.
#
#  Note: 'velocity' and 'constraint' may appear multiple times
#
# Output:
#   fi: the analysed field
#   err: error variance of the analysis field relative to the error variance of 
#     the background
#   s: structure
#     s.iB: adimensional
#     s.E: scaled EOF (dimensional)
#
# Note:
#   If zero is not a valid first guess for your variable (as it is the case for 
#   e.g. ocean temperature), you have to subtract the first guess from the 
#   observations before calling divand and then add the first guess back in.
#
# Example:
#   see divand_simple_example.m
#


function divandrun(mask,pmn,xi,x,f,len,lambda; 
                velocity = (),
                EOF = [],
                diagnostics = 0,
                EOF_lambda = 0,
                primal = 1,
                factorize = 1,
                tol = 1e-6,
                maxit = 100,
                minit = 10,
                constraints = (),
                inversion = "chol",
                moddim = [],
                fracindex = [],
                alpha = [],
                keepLanczosVectors = 0
#,
#                compPC = divand_pc_sqrtiB
                )


# check inputs

if !any(mask[:])
  error("no sea points in mask");
end

s = divand_background(mask,pmn,len,alpha,moddim);
s.betap = 0;
s.EOF_lambda = EOF_lambda;
s.primal = primal;
s.factorize = factorize;
s.tol = tol;
s.maxit = maxit;
s.minit = minit;
s.inversion = inversion;
s.keepLanczosVectors = keepLanczosVectors;
#s.compPC = compPC;

# remove non-finite elements from observations
f = f(:);
valid = isfinite(f);
x = cat_cell_array(x);

if !all(valid)  
  fprintf(1,"remove %d (out of %d) non-finite elements from observation vector\n",sum(!valid),numel(f));
  x = reshape(x,[length(f) s.n]);
  f = f[valid];
  x = reshape(x(repmat(valid,[1 s.n])),[length(f) s.n]);
  
  if !isempty(fracindex)
    fracindex = fracindex[:,valid];
  end

  if isscalar(lambda)
    # do nothing
  elseif isvector(lambda)
    lambda = lambda[valid];
  elseif ismatrix(lambda)
    lambda = lambda[valid,valid];
  end
end

apply_EOF_contraint = !(isempty(EOF) | all(EOF_lambda == 0));

s.mode = 1;

if !apply_EOF_contraint
    s.betap = 0;
else
    if s.mode==0
        s.betap = max(EOF_lambda)/s.coeff;  # units m^(-n)
    elseif s.mode==1
        s.betap = max(max(EOF_lambda)-1,0)/s.coeff;
    end
end

#assert(s.betap,0,1e-8)
# increase contraint on total enegery to ensure system is still positive defined
#s.betap
s.iB = s.iB + s.betap * s.WE'*s.WE;

# add observation constrain to cost function
s = divand_addc(s,divand_obs(s,xi,x,f,lambda,fracindex));

# add advection constraint to cost function
for i=1:length(velocity)
    #s = divand_advection(s,velocity);
    s = divand_addc(s,divand_constr_advec(s,velocity[i]));
end


# add all additional constrains
for i=1:length(constraints)
    s = divand_addc(s,constraints[i]);
end


if apply_EOF_contraint
    s = divand_eof_contraint(s,EOF_lambda,EOF);
end

# factorize a posteori error covariance matrix
# or compute preconditioner
s = divand_factorize(s);

if !apply_EOF_contraint
    fi,s = divand_solve(s,f);
else
    fi,s = divand_solve_eof(s,f);
end

return fi
# varargout{1} = fi;

# if nargout-diagnostics >= 2
#     err = divand_error(s);
#     varargout{2} = err;
# end

# if diagnostics
#     s.B = CovarIS(s.iB);
#     [s.Jb,s.Jo,s.Jeof,s.J] = divand_diagnose(s,fi,f);
#     s.valid = valid;
#     varargout{nargout} = s;
# end

end

# Copyright (C) 2008-2016 Alexander Barth <barth.alexander@gmail.com>
#
# This program is free software; you can redistribute it and/or modify it under
# the terms of the GNU General Public License as published by the Free Software
# Foundation; either version 2 of the License, or (at your option) any later
# version.
#
# This program is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
# FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
# details.
#
# You should have received a copy of the GNU General Public License along with
# this program; if not, see <http://www.gnu.org/licenses/>.

# LocalWords:  fi divand pmn len diag CovarParam vel ceil moddim fracdim


