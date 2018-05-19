"""


    qcvalues = divand_qc(fi,s,method);

# Input:

* `fi` : interpolated field from a `divandrun` execution
* `s`: corresponding structure returned by `divand`
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

If you cannot run `divandrun` but use `divandgo` (which does not provide a structure s at the output),
the latter provides `qcvalues` if you call `divandgo` with a keyword parameter `QCMETHOD=`

"""


function divand_qc(fi, s, method=0)

    info("Applying quality check based on the analysis")
    # For the moment, hardwired values
    # Make sure to work only with real observations
    switchvalue1=130;

    H = s.obsconstrain.H;
    R = s.obsconstrain.R;
    yo=s.yo;

    nd=size(s.obsout)[1];
    invlam=mean(diag(R));
    d0d=dot((1-s.obsout).*(s.yo),(s.yo));
    nrealdata=sum(1-s.obsout);

    meaneps2=(d0d/nrealdata) *invlam/(1+invlam);

    qcval=zeros(nd);

    residual=(1-s.obsout).*divand_residualobs(s,fi);



    if method==0

        mymethod=3
        if nrealdata < switchvalue1
            mymethod=1
        end
    else
        mymethod=method
    end


    # Third method
    if mymethod==4

        cvval=1
        qcval=residual.^2./(cvval*(diag(R)/invlam).*(1-divand_GCVKiiobs(s)).^2);
        return qcval
    end


    if mymethod==1

        qcval=residual.^2./(meaneps2*(diag(R)/invlam).*(1-divand_diagHKobs(s)));
        return qcval
    end

    if mymethod==3

        qcval=residual.^2./(meaneps2*(diag(R)/invlam).*(1-divand_GCVKiiobs(s)));
        return qcval
    end

	if mymethod==5

        qcval=residual.^2./(meaneps2*(diag(R)/invlam).*(1-divand_GCVKiiobs(s,-1;FIELD=fi)));
        return qcval
    end



    warn("divand_qc not defined for method  $method")

    return 0

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

# LocalWords:  fi divand pmn len diag CovarParam vel ceil moddim fracdim
