# A simple example of divand in 2 dimensions
# with observations from an analytical function.


# To test a single point
#Analysis in very large domain,
#diagnose value at distance 3*l

# Then the same but domain size just sligtly larger than 6*l with data in the center

# Then the same but with data at the border

using divand
using PyPlot

lpmsize=30
dsoverlsize=30
alpha=zeros(lpmsize,dsoverlsize)
varb=zeros(lpmsize,dsoverlsize)
varr=zeros(lpmsize,dsoverlsize)

lpm=collect(linspace(4,40,lpmsize))
lpmc=zeros(lpmsize)
dsoverl=collect(linspace(4,20,dsoverlsize))

for ii=1:lpmsize
    for jj=1:dsoverlsize

        len=1.0/dsoverl[jj]

        testpm=lpm[ii]/len

        isize=Int(ceil(testpm))

        @show isize

        xi=0
        mask=0
        epsilon2=0
        x=0
        f=0
        pm=0


        aj=zeros(300)
        vj=zeros(300)

        for j=1:300
            #for j=1:1
            alen=j/50

            alen=2

            # observations
            x = [0.5];
            f = [1];

            xi = collect(linspace(0,1,isize));



            # all points are valid points
            mask = trues(xi);

            # this problem has a simple cartesian metric
            # pm is the inverse of the resolution along the 1st dimension
            # pn is the inverse of the resolution along the 2nd dimension

            pm = ones(xi) / (xi[2]-xi[1]);
            # obs. error variance normalized by the background error variance
            epsilon2 = 10000;

            #Test to push boundary to wider distance:

            @show pm[1]*len
            @show len
            lpmc[ii]=pm[1]*len




            pm[isize]=1./(alen*len);
            pm[1]  =1./(alen*len);


            # correlation length



            # fi is the interpolated field
            fi2,s = divandrun(mask,(pm,),(xi,),(x,),f,len,epsilon2);


            #pcolor(reshape(diag(s.P),59,59)')
            #colorbar()

            aj[j]=alen
            vj[j]=var(diag(s.P))

        end
        alpha[ii,jj]=aj[indmin(vj)]
        varb[ii,jj]=vj[indmin(vj)]

        # now reference var
        pm = ones(xi) / (xi[2]-xi[1]);
        # obs. error variance normalized by the background error variance
        epsilon2 = 10000;

        #Test to push boundary to wider distance:

        @show pm[1]*len

        # fi is the interpolated field
        fi2,s = divandrun(mask,(pm,),(xi,),(x,),f,len,epsilon2);


        #pcolor(reshape(diag(s.P),59,59)')
        #colorbar()


        varr[ii,jj]=var(diag(s.P))



    end
end



pcolor(lpmc,dsoverl,varb)
colorbar()
clim(0,0.005)
savefig("var.withalpha2.png")

pcolor(lpmc,dsoverl,varr)
colorbar()
clim(0,0.1)
savefig("varr.png")


# Copyright (C) 2014, 2017 Alexander Barth <a.barth@ulg.ac.be>
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
