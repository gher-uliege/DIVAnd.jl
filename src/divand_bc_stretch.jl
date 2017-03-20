"""

"""
function divand_bc_stretch(mask,pmnin,xiin,Lin,moddim,alphabc=1)






    # number of dimensions
    n = ndims(mask)
    sz = size(mask)

    if isempty(moddim)
        moddim = zeros(1,n)
    end

    iscyclic = moddim .> 0

    Labs=deepcopy(Lin)
    #@show Labs
    if isa(Labs,Number)
        Labs = ((Labs * ones(size(mask)) for i=1:n)...)
    elseif isa(Labs,Tuple)

        if isa(Labs[1],Number)
            Labs = ([Labs[i] * ones(size(mask)) for i = 1:n]...)
        end

        for i=1:n
            if !isequal(size(mask),size(Labs[i]))
                error("mask $(formatsize(size(mask))) and correlation length $(formatsize(size(Labs[i]))) have incompatible size")
            end
        end
    end


    # Just used to fill the Labs tuple (so background will not fill again)
    #

    if alphabc==0
        warn("divand_bc_stretch was just used to fill in Labs")
        return pmnin,xiin,Labs
    end



    if alphabc>0

        pmn=deepcopy(pmnin)
        xi=deepcopy(xiin)

        for i=1:n



            #                 x first to work with undmodified pmninin


            if ~iscyclic[i]
                ind1 = [(j == i ? (1) : (:)) for j = 1:n]
                ind2 = [(j == i ? (2) : (:)) for j = 1:n]

                xi[i][ind1...]=xi[i][ind1...].+(xi[i][ind1...]-xi[i][ind2...]).*max( ((2.0.*alphabc.*Labs[i][ind1...].*pmnin[i][ind2...].-1.0).*pmnin[i][ind1...].-pmnin[i][ind2...])./(pmnin[i][ind1...]+pmnin[i][ind2...])    ,0.)
                ind1 = [(j == i ? (sz[i]) : (:)) for j = 1:n]
                ind2 = [(j == i ? (sz[i]-1) : (:)) for j = 1:n]

                xi[i][ind1...]=xi[i][ind1...]+(xi[i][ind1...]-xi[i][ind2...]).*max( ((2.0.*alphabc.*Labs[i][ind1...].*pmnin[i][ind2...]-1.0).*pmnin[i][ind1...]-pmnin[i][ind2...])./(pmnin[i][ind1...]+pmnin[i][ind2...])    ,0.)
                #                                       xi[i][end]=xi[i][end]+(xi[i][end]-xi[i][end-1])*max( ((2.0.*alphabc*Labs[i][end]*pmnin[i][end-1]-1.0)*pmnin[i][end]-pmnin[i][end-1])/(pmnin[i][end]+pmnin[i][end-1])    ,0.)
            end



            # if n==2
            # if i==1
            # if ~iscyclic[i]
            # #                                     @show max( ((2*alphabc*Labs[i][1,:].*pmnin[i][2,:]-1).*pmnin[i][1,:]-pmnin[i][2,:])./(pmnin[i][1,:]+pmnin[i][2,:])    ,0.0)
            # xi[i][1,:]=xi[i][1,:]+(xi[i][1,:]-xi[i][2,:]).*max( ((2*alphabc*Labs[i][1,:].*pmnin[i][2,:]-1).*pmnin[i][1,:]-pmnin[i][2,:])./(pmnin[i][1,:]+pmnin[i][2,:])    ,0.0)
            # xi[i][end,:]=xi[i][end,:]+(xi[i][end,:]-xi[i][end-1,:]).*max( ((2*alphabc*Labs[i][end,:].*pmnin[i][end-1,:]-1).*pmnin[i][end,:]-pmnin[i][end-1,:])./(pmnin[i][end,:]+pmnin[i][end-1,:])    ,0.0)
            # end
            # end
            # if i==2
            # if ~iscyclic[i]
            # xi[i][:,1]=xi[i][:,1]+(xi[i][:,1]-xi[i][:,2]).*max( ((2*alphabc*Labs[i][:,1].*pmnin[i][:,2]-1).*pmnin[i][:,1]-pmnin[i][:,2])./(pmnin[i][:,1]+pmnin[i][:,2])    ,0.0)
            # xi[i][:,end]=xi[i][:,end]+(xi[i][:,end]-xi[i][:,end-1]).*max( ((2*alphabc*Labs[i][:,end].*pmnin[i][:,end-1]-1).*pmnin[i][:,end]-pmnin[i][:,end-1])./(pmnin[i][:,end]+pmnin[i][:,end-1])    ,0.0)
            # end
            # end
            # end


            # now pmn



            wjmb=pmn[i]
            # For the moment, hardcoded for 1D and 2D

            #
            #                       if n==1
            if ~iscyclic[i]
                ind1 = [(j == i ? (1) : (:)) for j = 1:n]
                ind2 = [(j == i ? (2) : (:)) for j = 1:n]

                wjmb[ind1...]=1.0./max((2*alphabc.*Labs[i][ind1...].-1.0./wjmb[ind2...]),1.0./wjmb[ind2...])

                ind1 = [(j == i ? (sz[i]) : (:)) for j = 1:n]
                ind2 = [(j == i ? (sz[i]-1) : (:)) for j = 1:n]
                wjmb[ind1...]=1.0./max((2*alphabc.*Labs[i][ind1...].-1.0./wjmb[ind2...]),1.0./wjmb[ind2...])
            end
            #                       end

            # if n==2
            # if i==1
            # if ~iscyclic[1]
            # #                                         wjmb[1,:]=1.0./(alphabc.*Labs[1][1,:])
            # #                                 wjmb[end,:]=1.0./(alphabc.*Labs[1][end,:])
            # wjmb[1,:]=1.0./max((2*alphabc.*Labs[1][1,:].-1.0./wjmb[2,:]),1.0./wjmb[2,:])
            # wjmb[end,:]=1.0./max((2*alphabc.*Labs[1][end,:].-1.0./wjmb[end-1,:]),1.0./wjmb[end-1,:])
            # end
            # end
            # if i==2
            # if ~iscyclic[2]
            # #                                 wjmb[:,1]=1.0./(alphabc.*Labs[2][:,1])
            # #                             wjmb[:,end]=1.0./(alphabc.*Labs[2][:,end])
            # wjmb[:,1]=1.0./max((2*alphabc.*Labs[2][:,1].-1.0./wjmb[:,2]),1.0./wjmb[:,2])
            # wjmb[:,end]=1.0./max((2*alphabc.*Labs[2][:,end].-1.0./wjmb[:,end-1]),1.0./wjmb[:,end-1])
            # end
            # end
            # end


            #

        end

    end

    return pmn,xi,Labs
end

# Copyright (C) 2014, 2016 Alexander Barth <a.barth@ulg.ac.be>
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
