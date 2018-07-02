"""


    erri = DIVAnd_cpme_go(mask,pmn,xi,x,f,len,epsilon2; ...);



# Input:
*  Same arguments as DIVAndrun with in addition
*  `MEMTOFIT=`: keyword controlling how to cut the domain depending on the memory remaining available for inversion (not total memory)
*  `RTIMESONESCALES=` : if you provide a tuple of length scales, data are weighted differently depending on the numbers of neighbours they have. See `weight_RtimesOne` for details 


# Output:
*  `erri`: relative error field using the clever poor man's error approach. Result on the same grid as fi. `



 ONLY USE THIS VERSION IF YOU CANNOT RUN `DIVAndgo` with `:cmpe` activated (or directly `DIVAnd_cpme` if you can run `DIVAndrun`)


"""
function DIVAnd_cpme_go(mask,pmn,xi,x,f,Labs,epsilon2; otherargs...)


    
    errorscale=1;

    # The factor 1.70677 is the best one in 2D but should be slightly different for other dimensions
    # Could be a small improvement. Also used in DIVAnd_aexerr


    if isa(Labs,Tuple)
        len=([x./1.70766 for x in Labs]...,);
    else
        len=Labs./1.70766
    end


    
        cpme,bidon,zut =  DIVAndgo(mask,pmn,xi,x,ones(size(f)),len,epsilon2,:none; otherargs...);
		@show size(cpme),size(bidon),size(zut)
        cpme=errorscale .* max.(-cpme.+1,0)

    return cpme

end

# Copyright (C) 2008-2017 Alexander Barth <barth.alexander@gmail.com>
#                         Jean-Marie Beckers   <JM.Beckers@ulg.ac.be>
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

# LocalWords:  fi DIVAnd pmn len diag CovarParam vel ceil moddim fracdim
