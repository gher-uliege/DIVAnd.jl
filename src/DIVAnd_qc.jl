"""


    qcvalues = DIVAnd_qc(fi,s,method);

# Input:

* `fi` : interpolated field from a `DIVAndrun` execution
* `s`: corresponding structure returned by `DIVAnd`
* `method` : optional argument, which describes the method to be used:
 1  as for standard cross validation,
 3  as for GCV,
 4  with CV estimator to be used outside the routine,
 5  Poor man's GCV using data instead of random vector,
 0  automatic selection of method.


# Output

* `qcvalues`: quality check values, one for each data point.
The higher the value, the more suspect a data point is.
Absolute values of `qcvalues` might be not robust when analysis parameters are uncertain.
The ranking is however quite robust.

If you cannot run `DIVAndrun` but use `DIVAndgo` (which does not provide a structure s at the output),
the latter provides `qcvalues` if you call `DIVAndgo` with a keyword parameter `QCMETHOD=`

"""
function DIVAnd_qc(fi, s, method=0)

    # @info "Applying quality check based on the analysis"
    # For the moment, hardwired values
    # Make sure to work only with real observations
    switchvalue1=130;

    H = s.obsconstrain.H;
    R = s.obsconstrain.R;
    yo = s.yo;

    obsin = .!s.obsout

    nd = length(s.obsout)
    invlam = mean(diag(R)[obsin])

    d0d = s.yo[obsin] ⋅ s.yo[obsin]
    nrealdata = sum(obsin)
    meaneps2 = (d0d/nrealdata) * invlam/(1+invlam);

    @debug "meaneps2: meaneps2"

    qcval = zeros(nd);

    residual = DIVAnd_residualobs(s,fi);
    residual[s.obsout] .= 0

    if method==0
        mymethod=3
        if nrealdata < switchvalue1
            mymethod = 1
        end
    else
        mymethod = method
    end

    @debug "QC method: $mymethod"

    if mymethod == 1
        qcval .= residual.^2 ./ (meaneps2*(diag(R)/invlam).*(1 .- DIVAnd_diagHKobs(s)));
    elseif mymethod == 3
        qcval .= residual.^2 ./ (meaneps2*(diag(R)/invlam).*(1 .- DIVAnd_GCVKiiobs(s)));
    elseif mymethod == 4
        cvval = 1
        qcval .= residual.^2 ./ (cvval*(diag(R)/invlam).*(1 .- DIVAnd_GCVKiiobs(s)).^2);
    elseif mymethod == 5
        qcval .= residual.^2 ./ (meaneps2*(diag(R)/invlam).*(1 .- DIVAnd_GCVKiiobs(s,-1;FIELD=fi)));
    else
        @warn "DIVAnd_qc not defined for method  $method"
    end

    return qcval


end

# Copyright (C) 2008-2017 Alexander Barth <barth.alexander@gmail.com>
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
