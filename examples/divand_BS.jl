using divand
using PyPlot
using NetCDF

function loadbigfile(fname)

    data = readlines(open(fname,"r"))
    nobs = length(data)

    lon = zeros(nobs)
    lat = zeros(nobs)
    depth = zeros(nobs)
    time = Array{DateTime}(nobs)
    value = zeros(nobs)
    id = Array{String}(nobs)


    for i in 1:nobs
        rec = split(data[i])
        lon[i] = parse(Float64,rec[1])
        lat[i] = parse(Float64,rec[2])
        value[i] = parse(Float64,rec[3])
        depth[i] = parse(Float64,rec[4])
        time[i] = DateTime(rec[10])
        id[i] = rec[11]
    end

    return value,lon,lat,depth,time,id
end

include("../src/override_ssmult.jl")

#fname = joinpath(ENV["HOME"],"Data/Salinity.bigfile")
fname = joinpath("C:/JMB/BlackSea/","Salinity.bigfile")
#fname = "C:/JMB/BlackSea/Salinity.bigfile"
#bath_name = joinpath(ENV["HOME"],"Data/DivaData/Global/gebco_30sec_16.nc")
bath_name = joinpath("C:/JMB/BlackSea/","diva_bath.nc")
bath_name = joinpath("C:/JMB/BlackSea/","gebco_30sec_16.nc")
isglobal = true

if !isdefined(:value)
    value,lon,lat,depth,time,id = loadbigfile(fname)
end

@show size(value)

dx = dy = 0.1
dx = dy = 0.2
dx = dy = 0.02
lonr = 27:dx:42
latr = 40.4:dy:46.6

lonr = 27:dx:42
latr = 40:dy:47


depthr = [0., 10, 20, 30, 50, 75, 100, 125, 150, 200, 250, 300, 400, 500, 600, 700, 800, 900, 1000, 1100, 1200, 1300, 1400, 1500, 1750, 2000];
#depthr = [0, 10, 20, 30, 50, 75, 100];
#depthr = 0:10.:30.;


timer = 1:1.:12

mask,(pm,pn,po,pp),(xi,yi,zi,ti) = divand_rectdom(lonr,latr,depthr,timer)
@show size(mask)     ()

epsilon2 = 1

time2 = Dates.month.(time)

mxi,myi,mask2 = load_mask(bath_name,isglobal,minimum(lonr),maximum(lonr),dx,minimum(latr),maximum(latr),dy,depthr)

mask3 = repeat(mask2,inner = (1,1,1,length(timer)))


sz = size(mask)
@show sz

z = zeros(sz)
# correlation length in arc degree
lenx = fill(.3,sz)
leny = fill(.3,sz)
# correlation length in meters
lenz = 10 + zi/5
# correlation time-scale in month
lent = fill(1.,sz)


# Prepare background as mean vertical profile and time evolution. Just call divand in two dimensions forgetting x and y ...
#
 mask4=trues(size(mask3)[3],size(mask3)[4])
vm = mean(value)
va = value - vm
fm,sm=divandrun(mask4,(po[1,1,:,:],pp[1,1,:,:]),(zi[1,1,:,:],ti[1,1,:,:]),(depth,time2),va,(4*lenz[1,1,:,:],4*lent[1,1,:,:]),epsilon2)
vaanalyzed=sm.H*statevector_pack(sm.sv,(fm,))



fma=zeros(size(mask3))

for i=1:size(mask3)[1]
for ii=1:size(mask3)[2]
fma[i,ii,:,:]=fm[:,:]
end
end

n=ndims(mask3)

toaverage=[true true false false]

dimstokeep=[]
	for i=1:n
		if !toaverage[i]
		  dimstokeep=vcat(dimstokeep,[i])
		end
	end
    @show dimstokeep	
	

reshapeshape =([(toaverage[i] ? (1) : (size(mask3)[i])) for i = 1:n]...)
copyshape =([(toaverage[i] ? (size(mask3)[i]) : (1)) for i = 1:n]...)
@show reshapeshape
@show copyshape

fmaa=repeat(reshape(fm,reshapeshape), inner=copyshape)

#fmaa=repeat(reshape(fm,(1,1,size(mask3)[3],size(mask3)[4])), inner=(size(mask3)[1],size(mask3)[2],1,1))
@show size(fma)
va=va-vaanalyzed




@show var(fma-fmaa)


fmb,ffb=divand_averaged_bg(mask3,(pm,pn,po,pp),(xi,yi,zi,ti),(lon,lat,depth,time2),va,(lenx,leny,4*lenz,4*lent),epsilon2,toaverage)

@show var(fma-fmb)

pcolor(zi[1,1,:,:],ti[1,1,:,:],fm+vm)

#Todo: define new anomalies based on this profiles using sm.H of course. 
# Maybe an even better solution is to make a general averaging function



#@time fi,s = divandrun(mask3,(pm,pn,po,pp),(xi,yi,zi,ti),(lon,lat,depth,time2),va,(lenx,leny,z,lent),epsilon2)
#fi = fi + vm;

#@show mean(fi)

#@time fip,sp = divandrun(mask3,(pm,pn,po,pp),(xi,yi,zi,ti),(lon,lat,depth,time2),va,(1,1,0,0.00001),epsilon2)
#fip = fip+vm;

# tolerance on the gradient A x - b
tol = 1e-4
# tolerance on the result x
tolres = 1e-3

kwargs = [(:tol, tol),(:maxit,10000),(:minit,0)]


compPC(iB,H,R) = x -> s.P*x

# @time fi2,s = divandrun(mask3,(pm,pn,po,pp),(xi,yi,zi,ti),(lon,lat,depth,time2),va,(lenx,leny,lenz,lent),epsilon2;
#                        kwargs...,inversion=:pcg,operatortype=Val{:MatFun},fi0=fi
#                        ,compPC = compPC
#                        )
# fi2 = fi2+vm;

# @show s.niter


#divand_save(replace(@__FILE__,r".jl$",".nc"),mask,"salinity",fi)

@show mean(va)

@time figo,s = divandgo(mask3,(pm,pn,po,pp),(xi,yi,zi,ti),(lon,lat,depth,time2),va,(lenx,leny,z,lent),epsilon2)

@show mean(va)

fib = figo + vm+fma;

@show var(fib-fi)/var(fi)

divand_save(replace(@__FILE__,r".jl$","b.nc"),mask,"salinity",fib)

nothing