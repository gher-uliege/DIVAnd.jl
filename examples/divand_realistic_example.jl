using divand
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


fname = joinpath(ENV["HOME"],"Data/Salinity.bigfile")
bath_name = joinpath(ENV["HOME"],"Data/DivaData/Global/gebco_30sec_16.nc")
isglobal = true

if !isdefined(:value)
    value,lon,lat,depth,time,id = loadbigfile(fname)
end

dx = dy = 0.1
lonr = 27:dx:42
latr = 40:dy:47
depthr = 0:10.:30
timer = 1:1.:12

mask,(pm,pn,po,pp),(xi,yi,zi,ti) = divand_rectdom(lonr,latr,depthr,timer)

epsilon2 = 0.1

time2 = Dates.month.(time)

mxi,myi,mask2 = load_mask(bath_name,isglobal,minimum(lonr),maximum(lonr),dx,minimum(latr),maximum(latr),dy,depthr)

mask3 = repeat(mask2,inner = (1,1,1,length(timer)))


vm = mean(value)
va = value - vm
@time fi,s = divandrun(mask3,(pm,pn,po,pp),(xi,yi,zi,ti),(lon,lat,depth,time2),va,(1,1,0,0),epsilon2)
fi = fi+vm;


#@time fip,sp = divandrun(mask3,(pm,pn,po,pp),(xi,yi,zi,ti),(lon,lat,depth,time2),va,(1,1,0,0.00001),epsilon2)
#fip = fip+vm;

# tolerance on the gradient A x - b
tol = 1e-4
# tolerance on the result x
tolres = 1e-3

kwargs = [(:tol, tol),(:maxit,5000),(:minit,0)]


compPC(iB,H,R) = x -> s.P*x

@time fi2,s = divandrun(mask3,(pm,pn,po,pp),(xi,yi,zi,ti),(lon,lat,depth,time2),va,(1,1,10,0.12),epsilon2;
                        kwargs...,inversion=:pcg,operatortype=Val{:MatFun},fi0=fi
                        ,compPC = compPC
                        )
fi2 = fi2+vm;

@show s.niter

nothing
