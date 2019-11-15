"""
    fma,faanom = DIVAnd_averaged_bg(mask,pmn,xi,x,f,len,epsilon2,toaverage;moddim=[])

# Input:
As for DIVAndrun, including all dimensions before averaging

# additional argument:
* toaverage: Array of ndims of boolean telling if in the corresponding direction averaging must be done

# Presently NO optional arguments from DIVAndrun supported except moddim

# Output:

* fma: Analysis where in the directions where toaverage is true, the same value is found
* faanom: Data anomalies when the analysis is subtracted from the input field.

"""
function DIVAnd_averaged_bg(
    mask,
    pmn,
    xi,
    x,
    f,
    len,
    epsilon2,
    toaverage;
    moddim = [],
    filterbackground = 0,
)

    n = ndims(mask)

    if isempty(moddim)
        moddim = zeros(1, n)
    end


    if sum(toaverage) == n
        vm = mean(f)
        fma = fill(vm, size(mask))
        faanom = f .- vm
        return fma, faanom
    end

    if sum(toaverage) == 0
        @warn "no averaging was asked in averaging routine"
        fma, s = DIVAndrun(mask, pmn, xi, x, f, len, epsilon2)
        faanom = f - s.H * statevector_pack(s.sv, (fma,))
        return fma, faanom
    end






    # if average in a direction, just take the point 1 in this direction

    #    @show toaverage
    ind1 = [(toaverage[i] ? (1) : (:)) for i = 1:n]

    #       @show ind1

    #       @show trues(3)
    #       @show size(mask)



    # In the following maybe there are better ways to extract the relevant fields
    # instead of copying tuples and then extracting. Did look into things like
    # ww=( [ (toaverage[i] ? z[i]:()) for i=1:4]...)
    # but it left empty dimensions instead of taking them out
    #

    dimstokeep = []
    for i = 1:n
        if !toaverage[i]
            dimstokeep = vcat(dimstokeep, [i])
        end
    end
    #   @show dimstokeep

    xm = x[dimstokeep]
    moddimm = moddim[dimstokeep]

    pmnmf = pmn[dimstokeep]
    ximf = xi[dimstokeep]

    pmnm = ([x[ind1...] for x in pmnmf]...,)
    xim = ([x[ind1...] for x in ximf]...,)


    if isa(len, Number)
        lenm = len
    elseif isa(len, Tuple)
        if isa(len[1], Number)
            lenm = len[dimstokeep]
        else
            lenmf = len[dimstokeep]
            lenm = ([x[ind1...] for x in lenmf]...,)
        end
    end


    # For len need to check if nunber, tuple of numbers or tuples of tuples


    #maskmf=trues(mask)

    #@show typeof(maskmf)

    #maskm=maskmf[ind1]

    maskm = trues(size(xim[1]))
    #       @show size(maskm)

    #print("save DIVAndrun")
    #JLD2.@save "/tmp/DIVAndrun.jld2" maskm pmnm xim xm f lenm epsilon2 moddimm
    fm, sm = DIVAndrun(maskm, pmnm, xim, xm, f, lenm, epsilon2; moddim = moddimm, alphabc=0)

    fm = DIVAnd_filter3(fm, NaN, filterbackground)

    vaanalyzed = sm.H * statevector_pack(sm.sv, (fm,))
    vaanalyzed[sm.obsout] .= NaN
    faanom = f - vaanalyzed
    #@show extrema(vaanalyzed), sum(sm.obsout), extrema(faanom)

    reshapeshape = ntuple(i -> (toaverage[i] ? 1 : size(mask, i)), n)
    copyshape = ntuple(i -> (toaverage[i] ? size(mask, i) : 1), n)

    #       @show reshapeshape
    #       @show copyshape
    fma = repeat(reshape(fm, reshapeshape), inner = copyshape)


    return fma, faanom

end



# Copyright (C) 2008-2019 Alexander Barth <barth.alexander@gmail.com>
#                         Jean-Marie Beckers   <JM.Beckers@ulg.ac.be>
#
# This program is free software; you can redistribute it and/or modify it under
# the terms of the GNU General Public License as published by the Free Software
# Foundation; either version 2 of the License, or (at your option) any later
# version.
#
# This program is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
# FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
# details.
#
# You should have received a copy of the GNU General Public License along with
# this program; if not, see <http://www.gnu.org/licenses/>.

# LocalWords:  fi DIVAnd pmn len diag CovarParam vel ceil moddim fracdim
