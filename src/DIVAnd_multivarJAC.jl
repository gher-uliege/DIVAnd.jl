function DIVAnd_multivarJAC(mask,pmn,xi,x,f,lenin,epsilon2in;epsilon2jacobian=1.0,kwargs...)

"""

DIVAnd_multivarJAC(mask,pmn,xi,x,f,lenin,epsilon2in;epsilon2jacobian=1.0,kwargs...)
# Same input as DIVAnd_multivarEOF

# Additional keyword argument input:

* `epsilon2jacobian` : epsilon2 array (or value) which defines the strength of the multivariate constraint for each variable.

      If you use a very large value for one variable, it means that variable will not be modified by the coupling, but the other will be. So typically
      if you have a habitat variable which you want to influence a variable for which you have few data only. Then you assign a large epsilon2forjacobian to the habitat variable and
      a lower value to the other one.

# Output:

* `fi`: analysis

* `s` : structure

* `emap` : univariate error field

* `emapm` : multivariate error field

* `pseudovelocity` : the "advection" field that guided the analysis



#

"""



    if ndims(mask)<3
        @error "Need two dimensional fields and several variables"
    end

    epsilon2=deepcopy(epsilon2in)
    sz = size(mask)
    lxl=zeros(Float64,size(mask)[end])
    lyl=zeros(Float64,size(mask)[end])

    eps2jac=epsilon2jacobian
    if isa(eps2jac, Number)
        eps2jac=tuple(repeat([eps2jac],size(mask)[end])...)
    end
    @show eps2jac

    # Some saveguards here; len must be zero on last dimension
    # At the same time create the lenlow correlation length for the analysis
    # in the lower dimension (without "variable" dimension) for the eof amplitudes
    #
    len=deepcopy(lenin)

        variableL=false
        if isa(lenin, Number)
            @warn "Expecting  len array with zero in last direction; will be enforced"

            len=tuple(repeat([lenin],ndims(mask)-1)...,0.0)
            @debug "Len enforced $len"
        elseif isa(lenin, Tuple)

            if isa(lenin[1], Number)


               if lenin[end]!=0.0
                @warn "Expecting zero len in last direction; enforcing"
                #@show lenin[end]
                #@show lenin[1:end-1]
                len=tuple(lenin[1:end-1]...,0.0)
                @debug "len enforced, $len"
                end
            else
                variableL=true
                if sum(len[end])!=0.0
                @warn "Expecting zero len in last direction , enforcing"
                len=(len[1:end-1]...,0.0 .* len[1])
                end
            end

        end





    # univariate analysis
    fi,s=DIVAndrun(mask,pmn,xi,x,f,len,epsilon2;kwargs...)
    #




    emap,meth=DIVAnd_errormap(mask,pmn,xi,x,f,len,epsilon2,s;method=:cheap,kwargs...)
    sori=deepcopy(s)
    @show "error method in multivar $meth"

        vconstrain=[]
        ML=5
        pseudovelocity=[zeros(Float64,size(mask)) for idx in 1:ndims(mask)]
        for mloop=1:ML
            #
            # calculate pseudo-velocities. Maybe best to do it with already
            # existing operators
            # dfdy



            S = sparse_stagger(sz, 2, s.iscyclic[2])
            m = (S * mask[:]) .== 1
            (coucouc,) = statevector_unpack(s.sv, sparse_pack(mask) * S' * sparse_pack(m)' *s.Dx[2]*statevector_pack(s.sv, (fi,)))
            pseudovelocity[1] .= coucouc

            S = sparse_stagger(sz, 1, s.iscyclic[1])
            m = (S * mask[:]) .== 1
            (coucouc,) = statevector_unpack(s.sv, sparse_pack(mask) * S' * sparse_pack(m)' *s.Dx[1]*statevector_pack(s.sv, (fi,)))
            pseudovelocity[2].= -coucouc

            for jjj=3:ndims(mask)
               pseudovelocity[jjj].=zeros(Float64,size(mask))
            end
            # now scale each "2D layer" of the velocity ()

           uu=reshape(pseudovelocity[1],(prod(size(pseudovelocity[1])[1:2]),prod(size(pseudovelocity[1])[3:end])))
           vv=reshape(pseudovelocity[2],(prod(size(pseudovelocity[2])[1:2]),prod(size(pseudovelocity[2])[3:end])))

           for kk=1:size(uu)[2]
            vn=sqrt(sum(uu[:,kk].^2 .+ vv[:,kk].^2)/size(uu)[1])
            if vn>0
             uu[:,kk].=uu[:,kk]./vn
             vv[:,kk].=vv[:,kk]./vn
            end
           end



          uu=reshape(uu,(prod(size(pseudovelocity[1])[1:end-1]),prod(size(pseudovelocity[1])[end])))
          vv=reshape(vv,(prod(size(pseudovelocity[2])[1:end-1]),prod(size(pseudovelocity[2])[end])))
          lu=zeros(size(uu)[2])
          lv=zeros(size(uu)[2])
          wwww=ones(size(uu)[2])
           wtot=0
           for lll=1:size(uu)[2]
            wwww[lll]=sum(abs.(x[end].-lll).<0.4)
            wtot=wtot+wwww[lll]
           end
           wwww.=wwww./wtot
           #@show wwww


          #now for each physical point
           for ii=1:size(uu)[1]
            lu.=uu[ii,:]
            lv.=vv[ii,:]
            vvv=lu.^2 .+ lv.^2
            wwww.=1.0 ./max.(0.000001,emap[ii:size(uu)[1]:end])
            wwww.=wwww/sum(wwww)


            sorti=sortperm(vvv;rev=true)

            meanu=lu[sorti[1]]*wwww[sorti[1]]
            meanv=lv[sorti[1]]*wwww[sorti[1]]
            for kk=2:size(uu)[2]
                kks=sorti[kk]
                if meanu*lu[kks]+meanv*lv[kks]> 0
                    meanu=meanu+lu[kks]*wwww[kks]
                    meanv=meanv+lv[kks]*wwww[kks]

                else
                    meanu=meanu-lu[kks]*wwww[kks]
                    meanv=meanv-lv[kks]*wwww[kks]

                end
            end
            # apply possible weighting and scaling by L2 here
            if variableL
                lxl[:].= len[1][ii:size(uu)[1]:end]
                lyl[:].= len[2][ii:size(uu)[1]:end]
                else
                lxl[:].= len[1]
                lyl[:].= len[2]
            end



            uu[ii,:].=meanu.*sqrt.(lxl[:].*lyl[:]./eps2jac[:])
            vv[ii,:].=meanv.*sqrt.(lxl[:].*lyl[:]./eps2jac[:])

           end
            uu[isnan.(uu)].=0.0
            vv[isnan.(vv)].=0.0
            pseudovelocity[1]=reshape(uu,size(pseudovelocity[1]))
            pseudovelocity[2]=reshape(vv,size(pseudovelocity[2]))
            #
             vconstrain=DIVAnd_constr_advec(sori, tuple(pseudovelocity...))
            # pass the corresponding constraint to divandrun
             fi,s=DIVAndrun(mask,pmn,xi,x,f,len,epsilon2;constraints=[vconstrain],kwargs...)





        end
      # Error maps for the multivariate approach
      emapm,methm=DIVAnd_errormap(mask,pmn,xi,x,f,len,epsilon2,s;method=:cheap,constraints=[vconstrain],kwargs...)
      @show methm
    # Banzaii, finished
    return fi,s,emap,emapm,pseudovelocity
end