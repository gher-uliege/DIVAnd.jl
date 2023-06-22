"""
    error,method = DIVAnd_errormap(mask,pmn,xi,x,f,len,epsilon2,
    s;
    method = :auto,
    Bscale = false,
    otherargs...,);



# Input: same as for DIVAndrun WHICH MUST HAVE BEEN EXECUTED to get `s`

* `mask`: binary mask delimiting the domain. true is inside and false outside. For oceanographic application, this is the land-sea mask.

* `pmn`: scale factor of the grid. pmn is a tuple with n elements. Every
       element represents the scale factor of the corresponding dimension. Its
       inverse is the local resolution of the grid in a particular dimension.

*  `xi`: tuple with n elements. Every element represents a coordinate
  of the final grid on which the observations are interpolated

* `x`: tuple with n elements. Every element represents a coordinate of
  the observations

* `f`: value of the observations *minus* the background estimate (m-by-1 array).
    (see note)

* `len`: correlation length

* `epsilon2`: error variance of the observations (normalized by the error variance of the background field). `epsilon2` can be a scalar (all observations have the same error variance and their errors are decorrelated), a vector (all observations can have a difference error variance and their errors are decorrelated) or a matrix (all observations can have a difference error variance and their errors can be correlated). If `epsilon2` is a scalar, it is thus the *inverse of the signal-to-noise ratio*.

* `s`: this is the structure returned from the analysis itself.

# Optional input arguments specified as keyword arguments also as for DIVAnd

*`method` : the method to be used, valid are `:auto`, `:cheap`, `:precise`, `:cpme`, `:scpme`, `:exact`, `:aexerr`, `:diagapp`

            auto will select depenting on data coverage and length scale
            cheap will also select but restrict to cheaper methods
            precise will also select but prefer better approximations
            the other choices are just the ones among which the automatic choices will choose from. You can force the choice by specifying the method.

*`Bscale` : it `true` will try to take out the the boundary effects in the background error variance. Not possible with all methods

# Output:

* `error`: the error map

* `method`: the method used

"""



function DIVAnd_errormap(
    mask,
    pmn,
    xi,
    x,
    f,
    len,
    epsilon2,
    s;
    method = :auto,
    Bscale = false,
    rng=Random.GLOBAL_RNG,
    otherargs...
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


    if method == :auto


        # try to guess

        # small L
        # very low data coverage cpme
        # very high data coverage: scpme
        # otherwise: diagapp

        # larger L:
        # very high data coverage: scpme
        if smallL
            if Lowdata
                errmethod = :cpme
            else
                if Bigdata
                    errmethod = :scpme
                else
                    errmethod = :diagapp
                end
            end
        else
            if Bigdata
                errmethod = :scpme
            else
                errmethod = :aexerr
            end

        end

        # So best guess up to now
        # @show errmethod
    end



    if method == :cheap

        if smallL
            if Lowdata
                errmethod = :cpme
            else
                if Bigdata
                    errmethod = :scpme
                else
                    errmethod = :cpme
                end
            end
        else
            if Bigdata
                errmethod = :scpme
            else
                errmethod = :cpme
            end

        end
    end

    if method == :precise



        if smallL
            if Lowdata
                errmethod = :aexerr
            else
                if Bigdata
                    errmethod = :diagapp
                else
                    errmethod = :diagapp
                end
            end
        else
            if Bigdata
                errmethod = :diagapp
            else
                errmethod = :aexerr
            end
        end
    end


    if errmethod == :cpme && Bscale
        @warn "Sorry, that method does not allow rescaling by spatial dependance of B "
        ScalebyB = false
    end

    if errmethod == :scpme && Bscale
        @warn "Sorry, that method does not allow rescaling by spatial dependance of B "
        ScalebyB = false
    end

    if errmethod == :exact && Bscale
        # Or maybe if all info is there run locally ? Yes probably possible as aexerr also needs all infos ?
        @warn "You need to do that scaling by yourself, running diva again with a very high R matrix and divide by this second map"
        ScalebyB = false
    end

    if errmethod == :scpme && noP
        @warn "Sorry, that method needs s.P to be available. Will use cpme instead"
        errmethod = cpme
    end

    if errmethod == :exact && noP
        @warn "Sorry, that method needs s.P to be available. Will use aexerr instead"
        errmethod = aexerr
    end

    if errmethod == :diagapp && noP
        @warn "Sorry, that method needs s.P to be available. Will use aexerr instead"
        errmethod = aexerr
    end


    # Now calculate error depening on the method

    # @show errmethod, ScalebyB, pointsperbubblelimit

    if errmethod == :cpme
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

    if errmethod == :scpme
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
        DIVAnd_scalecpme!(scpme,s.P;rng=rng)

        return scpme, errmethod
    end

    if errmethod == :exact

        errormap, =statevector_unpack(s.sv,diag(s.P))

        return errormap, errmethod
    end

    if errmethod == :diagapp
        errormap = DIVAnd_diagapp(
            s.P,
            pmn,
            len,
            s.sv,
        )
        return errormap, errmethod
    end

    if errmethod == :aexerr
        errormap,bi,c,d =DIVAnd_aexerr(
            mask,
            pmn,
            xi,
            x,
            f,
            len,
            epsilon2;
            rng=rng,
            otherargs...
        )
        if errormap==0
        @warn "too fine resolution for aexerr, using exact"
        errormap, =statevector_unpack(s.sv,diag(s.P))

        return errormap, errmethod

        end
        if ScalebyB
        return errormap./bi, errmethod
        else
        return errormap, errmethod
        end
    end
    @show "You should not be here"
end
