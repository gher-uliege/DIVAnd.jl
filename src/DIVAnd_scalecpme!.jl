function DIVAnd_scalecpme!(cpme, P::CovarIS, nsamples = 7)
    # IN PLACE rescaling of the clever poor mans estimate using a randomized estimate of the analysis error covariance
    # P returned in structure s (so s.P) from a previous run
    # nsamples is the number of random arrays used to estimate the value

	fractionshift=0.5
    rng = Random.GLOBAL_RNG
	z = randn(rng,(size(P)[1], nsamples))
    errscale = 1
    if P.factors != nothing
        ZZ = P.factors.PtL \ z
        errscale = mean(diag(ZZ' * ZZ) ./ diag(z' * z))
    else
        ZZ = P * z
        errscale = mean(diag(z' * ZZ) ./ diag(z' * z))
    end
	# errscale here is the target

	# oldvalue is what we had
	oldvalue=mean(cpme[.!isnan.(cpme)])

	#try a shift of fraction fractionshift of the mismatch

	cpme[:] .= cpme[:] .+ fractionshift*(errscale-oldvalue)

	# the rest is corrected by the multiplicative scaling


    cpme[:] .= cpme[:] * errscale / mean(cpme[.!isnan.(cpme)])
    return cpme
end
