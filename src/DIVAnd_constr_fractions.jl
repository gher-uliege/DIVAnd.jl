"""
    c = DIVAnd_constr_fractions(s,epsfractions)

Creates integral constraints so that fractions sum up to "1". It is assumed that the different fractions to be analyzed are in the last dimension of the geometry

Input:
  s: structure
  
  epsfraction: error variance on constraint. 
  


Output:
  c: structure to be used by DIVAnd_addc with the following fields: R (a
    covariance matrix), H (extraction operator) and yo (specified value for
    the constrain).
"""
function DIVAnd_constr_fractions(s,epsfractions)

   

   
    

    sz = size(s.mask);
	
	nwithoutfractions=Int(floor((s.sv.n)/sz[end]))
	
	if mod(s.sv.n,sz[end])>0
	
	 warn("Problem in constraint for fractions, mask does seem to vary for the different fractions")
	
	end
	 




     A = spzeros(nwithoutfractions,s.sv.n)
	 l = size(A,1);
	 yo = ones(l)
	 R = Diagonal(epsfractions.*ones(l))


     for i=1:nwithoutfractions
	   for j=1:sz[end]
	   A[i,i+(j-1)*nwithoutfractions]=1
	   end
	 end
	   

    return DIVAnd_constrain(yo,R,A)

end

# Copyright (C) 2014, 2017 Alexander Barth <a.barth@ulg.ac.be>
#                     2018 Jean-Marie Beckers <JM.Beckers@uliege.be>
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
