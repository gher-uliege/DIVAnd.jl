function DIVAnd_errormap(
    mask,
    pmn,
    xi,
    x,
    f,
    len,
    epsilon2,
    s;
    method = "auto",
    Bscale = false,
    otherargs...,
)

    # Criteria to define which fraction of the domain size L can be to be called small
    LoverLlimit = 0.1
    # Criteria for lot of data means lot of data in hypersphere of correlation length diameters. Need to think about L=0 case ...
    pointsperbubblelimit = 10
    pointsperbubblelimitlow = 1


    errmethod = method
    ScalebyB = Bscale
    noP = s.P == ()

    smallL = false
    Bigdata = false
    Lowdata = false


    if method == "auto"

        Lpmnrange = DIVAnd_Lpmnrange(pmn, len)
        # L compared to domain size

        LoverLdomain = zeros(Float64, ndims(mask))


        for i = 1:ndims(mask)
            LoverLdomain[i] = Lpmnrange[i][2] / size(mask)[i]
        end

        if sum(LoverLdomain .< LoverLlimit) == ndims(mask)
            smallL = true
        end


        # Now look at lower values to check for data coverage
        realdims = ndims(mask)
        for i = 1:ndims(mask)
            LoverLdomain[i] = Lpmnrange[i][1] / size(mask)[i]


            if Lpmnrange[i][1] == 0
                LoverLdomain[i] = 1.0 / size(mask)[i]
                realdims = realdims - 1
            end
        end
        #nbdonnee size of f a revoir en fonction depsilon2
        if prod(LoverLdomain)*(pi^realdims)/gamma((realdims/2)+1) * size(f)[1] > pointsperbubblelimit
            Bigdata = true
        end

        if prod(LoverLdomain)*(pi^realdims)/gamma((realdims/2)+1) * size(f)[1] < pointsperbubblelimitlow
            Lowdata = true
        end

        # try to guess

        # small L
        # very low data coverage cpme
        # very high data coverage: scpme
        # otherwise: diagapp

        # larger L:
        # very high data coverage: scpme
        if smallL
            if Lowdata
                errmethod = "cpme"
            else
                if Bigdata
                    errmethod = "scpme"
                else
                    errmethod = "diagapp"
                end
            end
        else
		    if Bigdata
                errmethod = "scpme"
            else
                errmethod = "aexerr"
		    end

        end

        # So best guess up to now
        @show errmethod
    end



    if method == "cheap"
        Lpmnrange = DIVAnd_Lpmnrange(pmn, len)
        # L compared to domain size

        LoverLdomain = zeros(Float64, ndims(mask))


        for i = 1:ndims(mask)
            LoverLdomain[i] = Lpmnrange[i][2] / size(mask)[i]
        end

        if sum(LoverLdomain .< LoverLlimit) == ndims(mask)
            smallL = true
        end


        # Now look at lower values to check for data coverage
        realdims = ndims(mask)
        for i = 1:ndims(mask)
            LoverLdomain[i] = Lpmnrange[i][1] / size(mask)[i]
            if Lpmnrange[i][1] == 0
                LoverLdomain[i] = 1.0 / size(mask)[i]
                realdims = realdims - 1
            end
        end
        #nbdonnee size of f a revoir en fonction depsilon2
        if prod(LoverLdomain) * size(f)[1] > pointsperbubblelimit^realdims
            Bigdata = true
        end

        if prod(LoverLdomain) * size(f)[1] < pointsperbubblelimitlow^realdims
            Lowdata = true
        end
        if smallL
            if Lowdata
                errmethod = "cpme"
            else
                if Bigdata
                    errmethod = "scpme"
                else
                    errmethod = "cpme"
                end
            end
        else
		    if Bigdata
                errmethod = "scpme"
            else
                errmethod = "cpme"
		    end

        end
    end

    if method == "precise"
        Lpmnrange = DIVAnd_Lpmnrange(pmn, len)
        # L compared to domain size

        LoverLdomain = zeros(Float64, ndims(mask))


        for i = 1:ndims(mask)
            LoverLdomain[i] = Lpmnrange[i][2] / size(mask)[i]
        end

        if sum(LoverLdomain .< LoverLlimit) == ndims(mask)
            smallL = true
        end


        # Now look at lower values to check for data coverage
        realdims = ndims(mask)
        for i = 1:ndims(mask)
            LoverLdomain[i] = Lpmnrange[i][1] / size(mask)[i]
            if Lpmnrange[i][1] == 0
                LoverLdomain[i] = 1.0 / size(mask)[i]
                realdims = realdims - 1
            end
        end
        #nbdonnee size of f a revoir en fonction depsilon2
        if prod(LoverLdomain) * size(f)[1] > pointsperbubblelimit^realdims
            Bigdata = true
        end

        if prod(LoverLdomain) * size(f)[1] < pointsperbubblelimitlow^realdims
            Lowdata = true
        end


        if smallL
            if Lowdata
                errmethod = "aexerr"
            else
                if Bigdata
                    errmethod = "diagapp"
                else
                    errmethod = "diagapp"
                end
            end
        else
		    if Bigdata
                errmethod = "diagapp"
            else
                errmethod = "aexerr"
		    end
        end
    end


    if errmethod == "cpme" && Bscale
        warn("Sorry, that method does not allow rescaling by spatial dependance of B ")
        ScalebyB = false
    end

    if errmethod == "scpme" && Bscale
        warn("Sorry, that method does not allow rescaling by spatial dependance of B ")
        ScalebyB = false
    end

    if errmethod == "exact" && Bscale
        # Or maybe if all info is there run locally ? Yes probably possible as aexerr also needs all infos ?
        warn("You need to do that scaling by yourself, running diva again with a very high R matrix and divide by this second map")
        ScalebyB = false
    end

    if errmethod == "scpme" && noP
        warn("Sorry, that method needs s.P to be available. Will use cpme instead")
        errmethod = cpme
    end

    if errmethod == "exact" && noP
        warn("Sorry, that method needs s.P to be available. Will use aexerr instead")
        errmethod = aexerr
    end

    if errmethod == "diagapp" && noP
        warn("Sorry, that method needs s.P to be available. Will use aexerr instead")
        errmethod = aexerr
    end


    # Now calculate error depening on the method

    @show errmethod, ScalebyB, pointsperbubblelimit

    if errmethod == "cpme"
        errormap = DIVAnd_cpme(
            mask,
            pmn,
            xi,
            x,
            f,
            len,
            epsilon2;
            otherargs...
        )

        return errormap, errmethod
    end

    if errmethod == "scpme"
        errormap = DIVAnd_cpme(
            mask,
            pmn,
            xi,
            x,
            f,
            len,
            epsilon2;
            otherargs...
        )

        scpme=deepcopy(errormap)
        DIVAnd_scalecpme!(scpme,s.P)

        return scpme, errmethod
    end

    if errmethod == "exact"
        errormap =statevector_unpack(s.sv,diag(s.P))

        return errormap, errmethod
    end

    if errmethod == "diagapp"
        errormap = DIVAnd_diagapp(
            s.P,
            pmn,
            len,
            s.sv,
        )
        return errormap, errmethod
    end

    if errmethod == "aexerr"
        errormap,bi,c,d =DIVAnd_aexerr(
            mask,
            pmn,
            xi,
            x,
            f,
            len,
            epsilon2;
            otherargs...
        )
        return errormap, errmethod
    end
    @show "You should not be here"
end
