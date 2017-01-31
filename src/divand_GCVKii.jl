"""
Computes an estimate of the mean value of the diagonal of HK using GCV and the already solved analysisand it structure s

Kii = divand_GCVKii(s);

"""


function divand_GCVKii(s,nr=5)

#the second, optional argument is the number of random vectors nr used for the estimate


H = s.H;
R = s.R;



Z=randn(size(R)[1],nr);




   P = s.P;
   WW=P * (H'* (R \ Z));
   ZtHKZ=  Z'*(H*WW);
   ZtZ  =  Z'*Z;
# Now take average of the nr different estimates
   Kii=mean(diag(ZtHKZ)./diag(ZtZ));
return Kii

end

# Copyright (C) 2008-2017 Alexander Barth <barth.alexander@gmail.com>
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


