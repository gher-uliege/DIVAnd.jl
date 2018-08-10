"""
Computes the diagonal terms of the so called hat-matrix HK, using the already solved analysis and it structure s. Warning: might take some time

diagonalterms = DIVAnd_diagHK(s);

"""
function DIVAnd_diagHK(s)
    H = s.H;
    R = s.R;

    # to be replaced later exploiting the factors of P ?
    #
    P = s.P;
    #   WW=P * (H'* (R \ I));
    #   ItHKI =  I'*H*WW;

    # parenthesis to force the order of operations

    #ItHKI =  I' * (H * (P * (H' * (R \ I))));
    #    ItHKI =   (H * (P * (H' * (R \ I))));
    #    diagHKb = diag(ItHKI);

    diagHK =
        @static if VERSION >= v"0.7.0-beta.0"
            # workaround for bug
            # https://github.com/JuliaLang/julia/issues/27860
            diagLtCM(copy(H'),P,(H' * (R \ I)))
        else
            # copy should not be needed
            diagLtCM(H',P,(copy(H') * (R \ I)))
        end
    #    if (norm(diagHKb-diagHK)> norm(diagHK)*1E-7)
    #     @warn "WTF"
    #    end


    return diagHK

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

# LocalWords:  fi DIVAnd pmn len diag CovarParam vel ceil moddim fracdim
