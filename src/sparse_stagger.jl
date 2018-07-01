"""
    S = sparse_stagger(sz1,m,cyclic)

Create a sparse operator `S` for staggering a field in dimension `m`
for "collapsed" matrix of the size `sz1`.
`cyclic` is true if domain is cyclic along dimension m.
`false` is the default value
"""
sparse_stagger(sz1,m,cyclic = false) = _sparse_wsum(sz1,m,0.5,0.5,cyclic)

# Copyright (C) 2009,2018 Alexander Barth <a.barth@ulg.ac.be>
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; If not, see <http://www.gnu.org/licenses/>.
