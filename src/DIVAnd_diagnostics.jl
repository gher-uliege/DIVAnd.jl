"""
    norm1,norm2,norm3,epsilon2=DIVAnd_norms(fi,s)

### Input:

* `fi` and `s` from the corresponding `DIVANDrun` exectution  `fi,s = DIVAndrun(...`

### Output:

the different terms of the norm which is minimised

* `norm1`: the regularity terms

* `norm2`: the data-misfit terms

* `norm3`: the remaining constraints

* `epsilon2` : the average diagonal of `R`; that allows to look at the L-shape curve by looking at `Log(norm1)` vs `Log(norm2*epsilon2)`



"""
function DIVAnd_norms(fi,s)

    R = s.R
    iB = s.iB
    H = s.H
    fstate = statevector_pack(s.sv, (fi,))
    residuals=s.yo - (s.H) * fstate
    residualsobs=s.obsconstrain.yo - (s.obsconstrain.H) *fstate
    residualsobs[s.obsout] .= 0.0

    norm2= residualsobs'*(s.obsconstrain.R\residualsobs)
    norm3= residuals'*(R \ residuals)-norm2
    norm1=fstate'*(iB*fstate)
    dR=diag(s.obsconstrain.R)
    epsilon2=mean(dR[isfinite.(dR)])

    return norm1,norm2,norm3,epsilon2

end















