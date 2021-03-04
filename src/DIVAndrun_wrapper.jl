function DIVAndrun(
    mask,
    pmn,
    xyi,
    xsel,
    vaa,
    len_scaled,
    epsilon2,
    errortype;
    MEMTOFIT = nothing, # unused
    overlapfactor = nothing, # unused
    solver = :direct, # unused
    kwargs...,
)

    fi, s =
        @time DIVAnd.DIVAndrun(
            mask, pmn, xyi,
            xsel,
            vaa,
            len_scaled,
            epsilon2; alphabc = 0,
            kwargs...)

    qcdata = ()
    if haskey(Dict(pairs(kwargs)),:QCMETHOD)
        qc_method = kwargs.QCMETHOD
        qcdata = DIVAnd.DIVAnd_qc(fi, s, qc_method)
    end
    residual = DIVAnd.DIVAnd_residualobs(s, fi);
    # save memory
    s = nothing
    GC.gc()

    erri = ()
    if errortype == :cpme
        erri = DIVAnd.DIVAnd_cpme(mask,pmn,xyi,xsel,vaa,len_scaled,epsilon2; kwargs...)
    end


    scalefactore = DIVAnd.DIVAnd_adaptedeps2(vaa, residual, epsilon2, isnan.(residual))

    return fi, erri, residual, qcdata, scalefactore
end
