# Solve the variational problem.
# 
# fi = divand_solve(s,yo)
#
# Derive the analysis based on all contraints included in s and using the 
# observations yo
#
# Input:
#   s: structure created by divand_factorize
#   yo: value of the observations
#
# Output:
#   fi: analyzed field

function divand_solve!(s)

H = s.H;
sv = s.sv;
R = s.R;
# fix me: ignore argument
yo = s.yo;


#if s.primal
#    if strcmp(s.inversion,'chol')    
      P = s.P;
      fpi =  P * (H'* (R \ yo[:]));
#     else        
#       #H = double(H);
#       HiRyo = H'* (R \ yo(:));
      
#       fun = @(x) s.iB*x + H'*(R\ (H * x));
      
#       if ~s.keepLanczosVectors          
#           [fpi,s.flag,s.relres,s.iter] = pcg(fun,HiRyo,s.tol,s.maxit,s.M1,s.M2);
          
#           if s.flag ~= 0
#               warning('divand-noconvergence','pcg did not converge');
#           end
          
#       else
#           #pc = @(x) s.M2 \ (s.M1 \ x);
#           pc = [];
#           x0 = zeros(size(s.iB,1),1);
#           [fpi,Q,T,diag] = conjugategradient(fun,HiRyo,'tol',s.tol,...
#                                              'maxit',s.maxit,...
#                                              'minit',s.minit,...
#                                              'x0',x0,'renorm',1,'pc',pc);
#           s.iter = diag.iter;
#           s.relres = diag.relres;
#           s.P = CovarLanczos(Q,T);
#       end
#     end        
# else # dual
#     B = s.B;        
#     # fun(x) computes (H B H' + R)*x
#     fun = @(x) H * (B * (H'*x)) + R*x;    
    
#     [tmp,flag,relres,iter] = pcg(fun, yo,s.tol,s.maxit,s.funPC);
    
#     if (flag ~= 0)
#         error('divand:pcg', ['Preconditioned conjugate gradients method'...
#             ' did not converge %d %g %g'],flag,relres,iter);
#     end
    
#     fpi = B * (H'*tmp);
# end


fi, = statevector_unpack(sv,fpi);

fi[!s.mask] = NaN;

return fi
end

# Copyright (C) 2014,2017 Alexander Barth <a.barth@ulg.ac.be>
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
