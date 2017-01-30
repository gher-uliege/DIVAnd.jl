"""
Computes the cross validation estimator (d-hat(d))' inv(R) (d-hat(d)) / ( 1' inv(R) 1)
where the hat value is the analysis not using a data point

theta = divand_cvestimator(s,residual);

"""


function divand_cvestimator(s,residual)



return (residual'*(s.R\ residual))/ (ones(size(residual))'*(s.R\ ones(size(residual))))

end

# Copyright (C) 2008-2017 Alexander Barth <barth.alexander@gmail.com>
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


