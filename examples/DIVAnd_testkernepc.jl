#SBATCH --mem-per-cpu=16000


# A simple example of DIVAnd in 2 dimensions
# with observations from an analytical function.

using DIVAnd
using Compat: @info, range
using PyPlot
if VERSION >= v"0.7"
    using Statistics
end

include("./prep_dirs.jl")

x=[0.]
y=[0.]
z=[0.]
f=[1.]



mask,(pm,pn,po),(xi,yi,zi) = DIVAnd_rectdom(range(-1,stop=1,length=40),range(-2,stop=2,length=45),range(-1.5,stop=1.5,length=50))

# correlation length
len = (0.1,0.2,0.15)

# obs. error variance normalized by the background error variance
epsilon2 = 1.;

@time fi,s = DIVAndrun(mask,(pm,pn,po),(xi,yi,zi),(x,y,z),f,len,epsilon2;alphabc=2);

itest=20
jtest=20
rms=zeros(itest)
epsfac=zeros(itest)
rms2=zeros(itest)
epsfac2=zeros(itest)

ffac=zeros(itest)
ffac2=zeros(itest)
rmsb=zeros(itest)
rmsb2=zeros(itest)


lenf=zeros(itest)
fipc=0
fipc2=0

for ii=1:itest


    lenfac=0.5+ii/itest*4.0

    lenfac=4.1
    epsfacc=0.00005+(ii-1)/itest*0.02

    lenf[ii]=lenfac

    lena=([x[1]*lenfac for x in len]...,)
    @time fipc,s = DIVAndrun(mask,(pm,pn,po),(xi,yi,zi),(x,y,z),f,lena,epsilon2;alphabc=2,alpha=[1,2,1]);

    @time fipc2,s = DIVAndrun(mask,(pm,pn,po),(xi,yi,zi),(x,y,z),f,lena,epsilon2;alphabc=2,alpha=[1,1]);

    ampli=sum(fipc.*fi)/sum(fipc.*fipc)
    ampli2=sum(fipc2.*fi)/sum(fipc2.*fipc2)

    ffac[ii]=ampli
    ffac2[ii]=ampli2

    epsfac[ii]=epsilon2/(ampli*(1+epsilon2)-1)
    epsfac2[ii]=epsilon2/(ampli2*(1+epsilon2)-1)

    epsfac[ii]=epsfacc
    epsfac2[ii]=epsfacc


    rms[ii]=sqrt(var(fi-ampli*fipc))
    rms2[ii]=sqrt(var(fi-ampli2*fipc2))

    @time fipc,s = DIVAndrun(mask,(pm,pn,po),(xi,yi,zi),(x,y,z),f,lena,epsilon2*epsfacc;alphabc=2,alpha=[1,2,1]);

    @time fipc2,s = DIVAndrun(mask,(pm,pn,po),(xi,yi,zi),(x,y,z),f,lena,epsilon2*epsfacc;alphabc=2,alpha=[1,1]);

    rmsb[ii]=sqrt(var(fi-fipc))
    rmsb2[ii]=sqrt(var(fi-fipc2))


end



subplot(1,2,1)
plot(lenf,rms,"-",lenf,rmsb,".")

subplot(1,2,2)
plot(lenf,rms2,"-",lenf,rmsb2,".")

figname = joinpath(figdir,basename(replace(@__FILE__,r".jl$" => "_1.png")));
savefig(figname)
@info "Saved figure as " * figname

figure("next")


subplot(1,2,1)
plot(lenf,epsfac,"-")

subplot(1,2,2)
plot(lenf,epsfac2,"-")

figname = joinpath(figdir,basename(replace(@__FILE__,r".jl$" => "_2.png")));
savefig(figname)
@info "Saved figure as " * figname

figure("nextb")

subplot(1,2,1)
plot(lenf,ffac,"-")

subplot(1,2,2)
plot(lenf,ffac2,"-")

figname = joinpath(figdir,basename(replace(@__FILE__,r".jl$" => "_3.png")));
savefig(figname)
@info "Saved figure as " * figname

# Copyright (C) 2014, 2018 Alexander Barth <a.barth@ulg.ac.be>
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
