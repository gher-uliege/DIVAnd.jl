function ufill(c,valex)
    imax,jmax,kmax = size(c)
    work = zeros(imax+2, jmax+2, kmax+2)
    work2 = zeros(imax+2, jmax+2, kmax+2)

    ic = zeros(Int8,imax+2, jmax+2, kmax+2)
    ic2 = zeros(Int8,imax+2, jmax+2, kmax+2)
    #ic = falses(imax+2, jmax+2, kmax+2)
    #ic2 = falses(imax+2, jmax+2, kmax+2)

    cfilled = copy(c)
    ufill!(cfilled,valex,work,work2,ic,ic2)

    return cfilled
end

function ufill!(c,valexc,work,work2,ic::Array{Int8,3},ic2::Array{Int8,3})
    const A1 = 5
    const A2 = 0
    const A3 = 0

    # fill the borders...
    # write(6,*) ' Filling in',valexc,imax,jmax,kmax

    imax,jmax,kmax = size(c)

    for j = 1:jmax+2
        for i = 1:imax+2
            work[i,j,1] = valexc
            ic[i,j,1] = 0
            work[i,j,kmax+2] = valexc
            ic[i,j,kmax+2] = 0
        end
    end

    for k = 1:kmax+2
        for i = 1:imax+2
            work[i,1,k] = valexc
            ic[i,1,k] = 0
            work[i,jmax+2,k] = valexc
            ic[i,jmax+2,k] = 0
        end
    end

    for k = 1:kmax+2
        for j = 1:jmax+2
            work[1,j,k] = valexc
            ic[1,j,k] = 0
            work[imax+2,j,k] = valexc
            ic[imax+2,j,k] = 0
        end
    end

    #
    # copy interior field
    for k = 1:kmax
        for j = 1:jmax
            for i = 1:imax
                work[i+1,j+1,k+1] = c[i,j,k]
                ic[i+1,j+1,k+1] = 1
                if work[i+1,j+1,k+1] == valexc
                    ic[i+1,j+1,k+1] = 0
                end
            end
        end
    end

    icount  =  1
    isom1 = isom2 = isom3 = 0
    
    while icount > 0
        icount = 0

        for k = 2:kmax+1
            for j = 2:jmax+1
                for i = 2:imax+1

                    work2[i,j,k] = work[i,j,k]
                    ic2[i,j,k] = ic[i,j,k]

                    if ic[i,j,k] == 0
                        work2[i,j,k] = valexc
                        icount = icount+1
                        isom = 0
                        
                        if A1 != 0
                            isom1 = (
                                #           ic[i,j,k+1]+ic[i,j,k-1]
                                0.
                                +ic[i+1,j,k]+ic[i-1,j,k]
                                +ic[i,j+1,k]+ic[i,j-1,k])
                        end
                            
                        if A2 != 0
                            isom2 = (
                                ic[i+1,j+1,k+1]+ic[i+1,j+1,k-1]
                                +ic[i+1,j-1,k+1]+ic[i+1,j-1,k-1]
                                +ic[i-1,j+1,k+1]+ic[i-1,j+1,k-1]
                                +ic[i-1,j-1,k+1]+ic[i-1,j-1,k-1])
                        end
                        
                        if A3 != 0
                            isom3 = (
                                ic[i,j+1,k+1]+ic[i,j+1,k-1]
                                + ic[i,j-1,k+1]+ic[i,j-1,k-1]
                                + ic[i+1,j,k+1]+ic[i+1,j,k-1]
                                + ic[i-1,j,k+1]+ic[i-1,j,k-1]
                                + ic[i+1,j+1,k]+ic[i+1,j-1,k]
                                + ic[i-1,j+1,k]+ic[i-1,j-1,k])
                        end
                        
                        isom = isom1*A1+isom2*A2+isom3*A3

                        if isom != 0
                            # interpolate
                            rsom1 = (
                                #           ic[i,j,k+1]*work[i,j,k+1]
                                #          +ic[i,j,k-1]*work[i,j,k-1]
                                0.
                                +ic[i+1,j,k]*work[i+1,j,k]
                                +ic[i-1,j,k]*work[i-1,j,k]
                                +ic[i,j+1,k]*work[i,j+1,k]
                                +ic[i,j-1,k]*work[i,j-1,k])
                            rsom2 = (
                                ic[i+1,j+1,k+1]*work[i+1,j+1,k+1]
                                +ic[i+1,j+1,k-1]*work[i+1,j+1,k-1]
                                +ic[i+1,j-1,k+1]*work[i+1,j-1,k+1]
                                +ic[i+1,j-1,k-1]*work[i+1,j-1,k-1]
                                +ic[i-1,j+1,k+1]*work[i-1,j+1,k+1]
                                +ic[i-1,j+1,k-1]*work[i-1,j+1,k-1]
                                +ic[i-1,j-1,k+1]*work[i-1,j-1,k+1]
                                +ic[i-1,j-1,k-1]*work[i-1,j-1,k-1])
                            rsom3 = (
                                ic[i,j+1,k+1]*work[i,j+1,k+1]
                                +ic[i,j+1,k-1]*work[i,j+1,k-1]
                                +ic[i,j-1,k+1]*work[i,j-1,k+1]
                                +ic[i,j-1,k-1]*work[i,j-1,k-1]
                                +ic[i+1,j,k+1]*work[i+1,j,k+1]
                                +ic[i+1,j,k-1]*work[i+1,j,k-1]
                                +ic[i-1,j,k+1]*work[i-1,j,k+1]
                                +ic[i-1,j,k-1]*work[i-1,j,k-1]
                                +ic[i+1,j+1,k]*work[i+1,j+1,k]
                                +ic[i+1,j-1,k]*work[i+1,j-1,k]
                                +ic[i-1,j+1,k]*work[i-1,j+1,k]
                                +ic[i-1,j-1,k]*work[i-1,j-1,k])
                            work2[i,j,k] = (A1*rsom1+A2*rsom2+A3*rsom3)/isom
                            ic2[i,j,k] = 1
                        end
                    end
                end
            end
        end

        for k = 2:kmax+1
            for j = 2:jmax+1
                for i = 2:imax+1
                    work[i,j,k] = work2[i,j,k]
                    ic[i,j,k] = ic2[i,j,k]
                end
            end
        end
        #@show icount

    end

# copy interior points
for k = 1:kmax
    for j = 1:jmax
        for i = 1:imax
            c[i,j,k] = work[i+1,j+1,k+1]
        end
    end
end
end



#imax = 100
#jmax = 100

c = randn(10,10,20)
valex = -9999
c[3:5,6:10,1:4] = valex

cf = ufill(c,valex);

SSH = Dataset("/home/abarth/Utils/sossheig.nc")["sossheig"][:];
c = copy(SSH.data); c[SSH.na] = valex;
@time cf = ufill(c,valex);
#pcolor(cf[:,:,1]')
@show extrema(cf)

