# A simple example of divand in 2 dimensions
# with observations from an analytical function.

using divand
using PyPlot
srand(1234)
# observations
nobs=100
x = 0.01+0.98*rand(nobs);
y = 0.01+0.98*rand(nobs);
f = sin(x*6) .* cos(y*6);
#f=-1+2*x
# x=[0.5,0.75]
# y=[0.5,0.75]
# f=[1,1]

# final grid
xi,yi = ndgrid(linspace(0,1,950),linspace(0,1,830));

# reference field
fref = sin(xi*6) .* cos(yi*6);

# all points are valid points
mask = trues(xi);

# this problem has a simple cartesian metric
# pm is the inverse of the resolution along the 1st dimension
# pn is the inverse of the resolution along the 2nd dimension

pm = ones(xi) / (xi[2,1]-xi[1,1]);
pn = ones(xi) / (yi[1,2]-yi[1,1]);

# correlation length
len = 0.03;

# obs. error variance normalized by the background error variance
epsilon2 = 1;
vscale=0.001
vscale=0

# fi is the interpolated field
@time fiexOLD,s = divandrun(mask,(pm,pn),(xi,yi),(x,y),f,(len,0.5*len),epsilon2;velocity=(vscale*yi,-vscale*xi),alphabc=0);

@time fiOLD,s = divandgo(mask,(pm,pn),(xi,yi),(x,y),f,(len,0.5*len),epsilon2;velocity=(vscale*yi,-vscale*xi),alphabc=0);



figure("P1")

subplot(1,3,1)
title("Old version")
pcolor(xi,yi,fiOLD)
clim(-1,1)
subplot(1,3,2)
title("Reference")
pcolor(xi,yi,fiexOLD)
clim(-1,1)
subplot(1,3,3)
rms=sqrt(var(fiexOLD-fiOLD))
title("rms $rms")
pcolor(xi,yi,fiexOLD-fiOLD)
clim(-0.05,0.05)
plot(x,y,"k.");
colorbar()

@time fiex,s = divandrun(mask,(pm,pn),(xi,yi),(x,y),f,(len,0.5*len),epsilon2;velocity=(vscale*yi,-vscale*xi),alphabc=2.);

@time fi,s = divandgo(mask,(pm,pn),(xi,yi),(x,y),f,(len,0.5*len),epsilon2;velocity=(vscale*yi,-vscale*xi),alphabc=2.);

figure("Pp")

subplot(1,3,1)
title("NEW version")
pcolor(xi,yi,fi)
clim(-1,1)
subplot(1,3,2)
title("Reference")
pcolor(xi,yi,fiex)
clim(-1,1)
subplot(1,3,3)
rms=sqrt(var(fiex-fi))
title("rms $rms")
pcolor(xi,yi,fiex-fi)
clim(-0.05,0.05)
plot(x,y,"k.");
colorbar()








# Copyright (C) 2014, 2017 Alexander Barth         <a.barth@ulg.ac.be>
#                          Jean-Marie Beckers   <JM.Beckers@ulg.ac.be>
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
