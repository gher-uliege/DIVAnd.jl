"""

"""
function DIVAnd_bc_stretch(mask,pmnin,xiin,Lin,moddim,alphabc=1)

    # number of dimensions
    n = ndims(mask)
    sz = size(mask)

    if isempty(moddim)
        moddim = zeros(1,n)
    end

    iscyclic = moddim .> 0

    Labs = len_harmonize(Lin,mask)

    # Just used to fill the Labs tuple (so background will not fill it again)
    #

    if alphabc==0
        #@warn "DIVAnd_bc_stretch was just used to fill in Labs"
        return pmnin,xiin,Labs
    end

    if alphabc>0
        pmn=deepcopy(pmnin)
        xi=deepcopy(xiin)

        for i=1:n

            if ~iscyclic[i]
                ind1 = [(j == i ? (1) : (:)) for j = 1:n]
                ind2 = [(j == i ? (2) : (:)) for j = 1:n]

                xi[i][ind1...] = xi[i][ind1...].+(xi[i][ind1...]-xi[i][ind2...]).*max.( ((2.0.*alphabc.*Labs[i][ind1...].*pmnin[i][ind2...].-1.0).*pmnin[i][ind1...].-pmnin[i][ind2...])./(pmnin[i][ind1...]+pmnin[i][ind2...])    ,0.)
                
				
				ind1 = [(j == i ? (sz[i]) : (:)) for j = 1:n]
                ind2 = [(j == i ? (sz[i]-1) : (:)) for j = 1:n]

                xi[i][ind1...] = xi[i][ind1...] .+ (xi[i][ind1...]-xi[i][ind2...]) .* max.(
                    ((2.0.*alphabc.*Labs[i][ind1...].*pmnin[i][ind2...] .- 1.0).*pmnin[i][ind1...]-pmnin[i][ind2...])./(pmnin[i][ind1...]+pmnin[i][ind2...])    ,0.)
                
            end





            # UPDATE BY REFERENCE
            wjmb=pmn[i]

            if ~iscyclic[i]
                ind1 = [(j == i ? (1) : (:)) for j = 1:n]
                ind2 = [(j == i ? (2) : (:)) for j = 1:n]

                wjmb[ind1...]=1.0./max.((2*alphabc.*Labs[i][ind1...].-1.0./wjmb[ind2...]),1.0./wjmb[ind2...])

                ind1 = [(j == i ? (sz[i]) : (:)) for j = 1:n]
                ind2 = [(j == i ? (sz[i]-1) : (:)) for j = 1:n]
            
				wjmb[ind1...]=1.0./max.((2*alphabc.*Labs[i][ind1...].-1.0./wjmb[ind2...]),1.0./wjmb[ind2...])
            end


        end

    end

    return pmn,xi,Labs
end

# Copyright (C) 2014, 2017 Alexander Barth 		<a.barth@ulg.ac.be>
#                         Jean-Marie Beckers 	<jm.beckers@ulg.ac.be>
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
