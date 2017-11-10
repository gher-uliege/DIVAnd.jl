using Base.Test


coord = [0.0853508 0.939756; 0.784134 0.080227; 0.999551 0.784304; 0.636594 0.7699; 0.357327 0.891722; 0.101827 0.856188; 0.862349 0.0555934; 0.992086 0.97036; 0.702955 0.591252; 0.685006 0.23132]

#x = (coord[:,1],coord[:,2])



LS = (0.1,0.2)

dist2(x,y,len) = sum(((x-y)./len).^2)
nobs = size(coord,1) 
x = ones(nobs)
Rx1 = zeros(nobs) 
Rx = zeros(nobs) 

function Rtimesx1!(coord,LS,x,Rx)
    nobs = size(coord,1) 
    len = [LS...]
    cov = zeros(nobs,nobs)
    
    for j = 1:nobs
        for i = 1:nobs
            d2 = dist2(coord[i,:],coord[j,:],len)
            cov[i,j] = exp(-d2)
        end
    end
    
    Rx[:] = cov*x
end


function Rtimesx!(coord,LS,x,Rx)
    nobs = size(coord,1) 
    len = [LS...]
    cov = zeros(nobs,nobs)

    Rx[:] = 0
    for j = 1:nobs
        for i = 1:nobs
            d2 = dist2(coord[i,:],coord[j,:],len)
            cov[i,j] = exp(-d2)

            Rx[i] += cov[i,j]*x[j]
        end
    end
    
end


@time Rtimesx1!(coord,LS,x,Rx1)
@time Rtimesx!(coord,LS,x,Rx)


@test Rx1 ≈ Rx
