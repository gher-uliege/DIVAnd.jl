"""
    ulon,ulat,meanval,stdval,count = statpos(val,lon,lat)

Return unique positions (`ulon`, `ulat`) as well their mean, 
standard deviation and count of the vector of observations `val` located
at the positions `lon` and `lat`.

"""

function statpos(val,lon,lat)
    allpos = collect(zip(lon,lat))
    uniquepos = collect(Set(allpos))

    count = zeros(Int,length(uniquepos))
    sumval = zeros(length(uniquepos))
    sumval2 = zeros(length(uniquepos))

    for i = 1:length(uniquepos)
        for j = 1:length(allpos)
            if uniquepos[i] == allpos[j]
                count[i] += 1
                sumval[i] += val[j]
                sumval2[i] += val[j]^2
            end
        end
    end
    meanval = sumval./count
    stdval = sumval2./count - meanval.^2;
    stdval[stdval .< 0] = 0;
    stdval = sqrt.(stdval)

    ulon = [p[1] for p in uniquepos];
    ulat = [p[2] for p in uniquepos];

    return ulon,ulat,meanval,stdval,count
end
