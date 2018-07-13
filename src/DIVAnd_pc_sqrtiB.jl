"""
Compute a preconditioner using the Cholesky decomposition.

[M1,M2] = DIVAnd_pc_michol(iB,H,R)

Compute preconditioner matrices M1 and M2 based on
the Cholesky decomposition of iB. The matrices H and R are not used.
M2 is the transpose of M1 for this preconditioner.
"""
function DIVAnd_pc_sqrtiB(iB,H,R)
    F =
        @static if VERSION >= v"0.7.0-beta.0"
            cholesky(iB);
        else
            cholfact(iB);
        end

    function fun!(x,fx)
        fx[:] = F \ x
    end
    return fun!
end

# LocalWords:  preconditioner diavnd pc michol iB Cholesky chol DIVAnd sqrtiB

# Copyright (C) 2017 Alexander Barth <a.barth@ulg.ac.be>
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
