function DIVAnd_scalecpme!(cpme,P::CovarIS,nsamples=7)
    # IN PLACE scaling of the clever poor mans estimate using a randomized estimate of the analysis error covariance
    # P returned in structure s (so s.P) from a previous run
	# nsamples is the number of random arrays used to estimate the value
    z=randn((size(P)[1],nsamples))
    errscale=1
    if P.factors != nothing
     ZZ=P.factors.PtL \ z
     errscale=mean(diag(ZZ'*ZZ)./diag(z'*z))
    else
     ZZ=P*z   
     errscale=mean(diag(z'*ZZ)./diag(z'*z))   
    end
    errscale=errscale/mean(cpme[.!isnan.(cpme)])
    cpme[:].=cpme[:]*errscale
    return cpme
end