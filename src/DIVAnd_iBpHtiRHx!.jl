function DIVAnd_iBpHtiRHx!(s,x::Array{Float64,1},iBx::Array{Float64,1},workobs1::Array{Float64,1},workstate1::Array{Float64,1}
      ,workstate2::Array{Float64,1},iBx_::Array{Float64,1};btrunc=Int[])


# Initialize outside workobs1, workstate1, workstate2,iBx_


	#iBx=s.iB*x
	mul!(iBx::Array{Float64,1},s.iB::SparseMatrixCSC{Float64,Int},x::Array{Float64,1})

    if btrunc==[]
		#iBx[:]=iBx[:]+s.H'*(s.R \ (s.H * x))
		mul!(workobs1::Array{Float64,1},s.H::SparseMatrixCSC{Float64,Int},x::Array{Float64,1})
		workobs1[:]=s.R\workobs1
        @static if VERSION >= v"0.7.0-beta.0"
		    mul!(workstate1::Array{Float64,1},(s.H::SparseMatrixCSC{Float64,Int})',workobs1::Array{Float64,1})
        else
		    At_mul_B!(workstate1::Array{Float64,1},s.H::SparseMatrixCSC{Float64,Int},workobs1::Array{Float64,1})
        end
		iBx[:]=iBx[:]+workstate1[:]
		return iBx
    end



    coeff = s.coeff;
    n = s.n;
    alpha=s.alpha

    # incomplete Bi calculated before now complemented on the fly


    for j=btrunc+1:length(alpha)
        k = Int(floor((j-2)/2))



        if mod(j,2) == 0


			#  D^k*x
			workstate1[:]=x[:]
			for kk=1:k
				mul!(workstate2::Array{Float64,1},s.D::SparseMatrixCSC{Float64,Int},workstate1::Array{Float64,1})
				workstate1[:]=workstate2[:]
			end


			# JMB Dirty Hack, stored Sum of (s.Dx[i]'*(s.WEss[i] *(s.WEss[i] *(s.Dx[i])))) into s.WEss[1] before

			    #                Dx = s.WEss[i] * (s.Dx[i] * Dk);
                #                iBx_ = iBx_ + Dx'*(Dx*x);

			mul!(iBx_::Array{Float64,1},s.WEss[1]::SparseMatrixCSC{Float64,Int},workstate1::Array{Float64,1})


			for kk=1:k
				#iBx_=Dk'*iBx_
                @static if VERSION >= v"0.7.0-beta.0"
				    mul!(workstate2::Array{Float64,1},(s.D::SparseMatrixCSC{Float64,Int})',iBx_::Array{Float64,1})
                else
				    At_mul_B!(workstate2::Array{Float64,1},s.D::SparseMatrixCSC{Float64,Int},iBx_::Array{Float64,1})
                end
				iBx_[:]=workstate2[:]
			end


        else


           # WD = s.WE * s.D^(k+1);


		   workstate1[:]=x[:]
		   for kk=1:k+1
		     mul!(workstate2::Array{Float64,1},s.D::SparseMatrixCSC{Float64,Int},workstate1::Array{Float64,1})
		     workstate1[:]=workstate2[:]
		   end

		   #iBx_ = WD'*(WD*x);

				mul!(workstate2::Array{Float64,1},s.WE::SparseMatrixCSC{Float64,Int},workstate1::Array{Float64,1})
				mul!(workstate1::Array{Float64,1},s.WE::SparseMatrixCSC{Float64,Int},workstate2::Array{Float64,1})

				# or better ?
				# workstate1[:]=(diag(s.WE).^2).*workstate1

			for kk=1:k+1
                @static if VERSION >= v"0.7.0-beta.0"
		            mul!(workstate2::Array{Float64,1},(s.D::SparseMatrixCSC{Float64,Int})',workstate1::Array{Float64,1})
                else
		            At_mul_B!(workstate2::Array{Float64,1},s.D::SparseMatrixCSC{Float64,Int},workstate1::Array{Float64,1})
                end
		     workstate1[:]=workstate2[:]
		    end
			iBx_[:]=workstate1[:]
        end



        #iBx_ = iBx_/coeff;
		#iBx = iBx + alpha[j] * iBx_

		asurc=alpha[j]/coeff
        iBx=BLAS.axpy!(asurc,iBx_,iBx)

    end


	#iBx[:]=iBx[:]+s.H'*(s.R \ (s.H * x))
	mul!(workobs1::Array{Float64,1},s.H::SparseMatrixCSC{Float64,Int},x::Array{Float64,1})
	workobs1[:]=s.R\workobs1
    if VERSION >= v"0.7.0-beta.0"
	    mul!(workstate1::Array{Float64,1},(s.H::SparseMatrixCSC{Float64,Int})',workobs1::Array{Float64,1})
    else
	    At_mul_B!(workstate1::Array{Float64,1},s.H::SparseMatrixCSC{Float64,Int},workobs1::Array{Float64,1})
    end
	iBx[:]=iBx[:]+workstate1[:]


    return iBx
end

# LocalWords:  iB DIVAnd

# Copyright (C) 2014-2017 Alexander Barth	  	 <a.barth@ulg.ac.be>
#                         Jean-Marie Beckers 	 <JM.Beckers@ulg.ac.be>
#
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
