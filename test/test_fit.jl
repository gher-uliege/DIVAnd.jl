function fitcorr(x,v,distbin,min_count)
    maxpoints = 100000

    pmax = length(distbin)-1;
    covar = zeros(pmax);
    count = zeros(pmax);
    distx = (distbin[1:end-1] + distbin[2:end])/2;
    
    distfun(xi,xj) = sqrt(sum(abs2,xi-xj))
    
    for l=1:maxpoints
    
        # random index
        i = rand(1:length(x[1]));
        j = rand(1:length(x[1]));

        xi = [c[i] for c in x]
        xj = [c[j] for c in x]

        if isnan(v[i]) || isnan(v[j])
            # one point is masked
            continue    
        end

        #dist = distance(y(i1,j1),x(i1,j1),y(i2,j2),x(i2,j2));
        distance = distfun(xi,xj)
  
        if distance >= distbin[end]
            # distance too large
            continue    
        end
  
        p = findlast(distance .>= distbin)
    
        if count[p] >= min_count
            # already enought points
            continue    
        end
  
        #distbin(p) <= dist && dist < distbin(p+1)
    
        covar[p] = covar[p] + v[i] * v[j]
        count[p] = count[p] + 1;

        #if mod(l,1000) == 0
        #    @show count
        #end
  
        if all(count .>= min_count)
            break
        end
    end
    

    covar = covar./count;

    return distx,covar
end

min_count = 5000

@time distx,covar = fitcorr(x,v,distbin,min_count)
