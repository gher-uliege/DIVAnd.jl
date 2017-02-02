"""
Compute a variational analysis of arbitrarily located observations to calculate data quality estimators

qcvalues,indexes = divand_qc(fi,s);


If the optional argument cvval is used, it should be the optimal value provided by cross validation

"""


function divand_qc(fi, s, cvval=0)


# For the moment, hardwired values
# Make sure to work only with real observations
switchvalue=500;

H = s.obsconstrain.H;
R = s.obsconstrain.R;
yo=s.yo;

nd=size(s.obsout)[1];
invlam=mean(diag(R));

meaneps2=(dot(yo[1:nd]',yo[1:nd])/nd) *invlam/(1+invlam);

qcval=zeros(nd);

residual=(1-s.obsout).*divand_residualobs(s,fi);
nrealdata=sum(1-s.obsout);

# Third method
if cvval>0

   c1=3
   qcval=residual.^2./(cvval*(diag(R)/invlam).*(1-divand_GCVKiiobs(s)).^2);
   return qcval,meaneps2,c1,invlam
end


if nrealdata<switchvalue

    c1=1
#   cvval=divand_cvestimator(s,residual./(1-divand_diagHKobs(s)));

   qcval=residual.^2./(meaneps2*(diag(R)/invlam).*(1-divand_diagHKobs(s)));
   else

    c1=2
#   cvval=divand_cvestimator(s,residual./(1-divand_GCVKiiobs(s)));	 
   qcval=residual.^2./(meaneps2*(diag(R)/invlam).*(1-divand_GCVKiiobs(s)));
end



return qcval,meaneps2,c1,invlam

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


