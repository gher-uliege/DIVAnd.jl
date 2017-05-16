"""



"""

function divand_iBx!(s,x::Array{Float64,1},iBx::Array{Float64,1},workobs1::Array{Float64,1},workstate1::Array{Float64,1}
      ,workstate2::Array{Float64,1},iBx_::Array{Float64,1},workdx::Array{Float64,1};btrunc=Int[])


# Initialize outside workobs1, workstate1, workstate2,iBx_
    #gc()
    #iBx[:]=s.iB*x   ::Array{Float64,1}
    #@show typeof(s.iB)
	A_mul_B!(iBx::Array{Float64,1},s.iB::SparseMatrixCSC{Float64,Int64},x::Array{Float64,1})
	
    if btrunc==[]
	    A_mul_B!(workobs1::Array{Float64,1},s.H::SparseMatrixCSC{Float64,Int64},x::Array{Float64,1})
		workobs1[:]=s.R\workobs1
		At_mul_B!(workstate1::Array{Float64,1},s.H::SparseMatrixCSC{Float64,Int64},workobs1::Array{Float64,1})
		iBx[:]=iBx[:]+workstate1[:]
		
	    #iBx[:]=iBx[:]+s.H'*(s.R \ (s.H * x))
        return iBx
    end

    
    #WE = s.WE;
    coeff = s.coeff;
    n = s.n;
    alpha=s.alpha
    # incomplete Bi calculated before now complemented on the fly

    #D=s.D
    for j=btrunc+1:length(alpha)
        # exponent of laplacian
        k = Int(floor((j-2)/2))
        # @show size(iBx)
        # @show size(x)
        # @show size(s.iB)
        #iBx_ = spzeros(iBx);
        iBx_[:]=0.*iBx ::Array{Float64,1};
# Certainly a gain to make; not recompute D^(k+1) but Dk*D if k+1 is one larger than already calculated value if it exists
# But only for n larger than 5 probably, so not an urgent thing
        
        if mod(j,2) == 0
            # constrain of derivative with uneven order (j-1)
            # (gradient, gradient*laplacian,...)
            # normalized by surface
			
			
			# Make a loop instead of matrix power !
			
			workstate1[:]=x[:]
			for kk=1:k
				A_mul_B!(workstate2::Array{Float64,1},s.D::SparseMatrixCSC{Float64,Int64},workstate1::Array{Float64,1})
				workstate1[:]=workstate2[:]
			end
			# contains D^k*x
			
			
			
			
			#if k>0
			#@show size(s.D)
			#Dk=(s.D^k) #::SparseMatrixCSC{Float64,Int64}
			#Dkx=Dk*x ::Array{Float64,1}
			#    else
			#Dkx=x ::Array{Float64,1}
			#end
			
            #for i=1:n
			# Dirty Hack 
			for i=1:1
			#? Not take up ?
			    if s.Ld[i]>0
                #                Dx = s.WEss[i] * (s.Dx[i] * Dk);
                #                iBx_ = iBx_ + Dx'*(Dx*x);
				#@show i
                #@time Dx = s.WEss[i] * s.Dx[i];
                # maybe gain if Dk=Dk'? Check with alex if s.Wess is diagonal
                #iBx_ = iBx_ + Dk'*(Dx'*(Dx*Dkx));
				#iBx_ = iBx_ + (s.Dx[i]'*(s.WEss[i] *(s.WEss[i] *(s.Dx[i]*Dkx))));
				
				#
				#@show size(workstate2), size(s.Dx[i]), size(workstate1)
				#A_mul_B!(workdx::Array{Float64,1},s.Dx[i]::SparseMatrixCSC{Float64,Int64},workstate1::Array{Float64,1})
				# if s.Dx[i] have different sizes, use local array
				
				#wwdx=s.Dx[i]*workstate1
				
				# s.WEss[i] is a diagonal matrix so in reality one just needs to scale workstate2
				
				
				#wwdx[:]=diag(s.WEss[i]).*(diag(s.WEss[i]).*wwdx[:])
				
				
				
				#workdx[:]=diag(s.WEss[i]).*workdx[:]
				#workdx[:]=diag(s.WEss[i]).*workdx[:]
				
				#A_mul_B!(workstate1::Array{Float64,1},s.WEss[i]::SparseMatrixCSC{Float64,Int64},workstate2::Array{Float64,1})
				
				
				# With dirtyhack
				A_mul_B!(workstate2::Array{Float64,1},s.WEss[i]::SparseMatrixCSC{Float64,Int64},workstate1::Array{Float64,1})
				
				#At_mul_B!(workstate1::Array{Float64,1},s.Dx[i]::SparseMatrixCSC{Float64,Int64},workdx::Array{Float64,1})
				
				# Try to see
				#At_mul_B!(workstate1::Array{Float64,1},s.Dx[i]::SparseMatrixCSC{Float64,Int64},wwdx::Array{Float64,1})
				
				iBx_[:] = iBx_[:] + workstate2[:]
				
				#iBx_ = iBx_ + (s.Dx[i]'*(s.WEss[i] *(s.WEss[i] *(s.Dx[i]*Dkx))));
				
				
				#iBx_ = iBx_ + (Dx'*(Dx*Dkx));
				end
			end
			
			for kk=1:k
			    
			    
				#iBx_=Dk'*iBx_
				At_mul_B!(workstate2::Array{Float64,1},s.D::SparseMatrixCSC{Float64,Int64},iBx_::Array{Float64,1})
				iBx_[:]=workstate2[:]
			end
        else
            # constrain of derivative with even order (j-1)
            # (laplacian, biharmonic,...)

            # normalize by surface of each cell
            # such that inner produces (i.e. WE'*WE)
            # become integrals
            # WD: units length^(n/2)
           #@show k 
           #@time WD = s.WE * s.D^(k+1);
		   
		   
		   workstate1[:]=x[:]
		   for kk=1:k+1
		     A_mul_B!(workstate2::Array{Float64,1},s.D::SparseMatrixCSC{Float64,Int64},workstate1::Array{Float64,1})
		     workstate1[:]=workstate2[:]
		   end
		   
		   #if k>-1
		   #@show k, size(s.D),typeof(s.D),s.D[1,1],typeof(s.D^(k+1))
		   # sDkp=(s.D^(k+1)) #::SparseMatrixCSC{Float64,Int64}
			#iBx_ = WD'*(WD*x);
			
				A_mul_B!(workstate2::Array{Float64,1},s.WE::SparseMatrixCSC{Float64,Int64},workstate1::Array{Float64,1})
				A_mul_B!(workstate1::Array{Float64,1},s.WE::SparseMatrixCSC{Float64,Int64},workstate2::Array{Float64,1})
				
				# or better ?
				# workstate1[:]=(diag(s.WE).^2).*workstate1
				
			for kk=1:k+1
		     At_mul_B!(workstate2::Array{Float64,1},s.D::SparseMatrixCSC{Float64,Int64},workstate1::Array{Float64,1})
		     workstate1[:]=workstate2[:]
		    end	
			iBx_[:]=workstate1[:]
				
			
			
			
		#	iBx_ = sDkp'*(s.WE*(s.WE*(sDkp*x)))
		 #	   else
		#	iBx_ = (s.WE*(s.WE*(x)))   
		#   end
			
			#Dk=s.WE * s.D^(k+1)
            #Dkx=Dk*x
            #iBx_ = Dk'*Dkx
        end
        asurc=alpha[j]/coeff
        #iBx_ = iBx_/coeff;
		#iBx = iBx + alpha[j] * iBx_

        iBx=BLAS.axpy!(asurc,iBx_,iBx)
        
    end
	
	A_mul_B!(workobs1::Array{Float64,1},s.H::SparseMatrixCSC{Float64,Int64},x::Array{Float64,1})
	workobs1[:]=s.R\workobs1
	At_mul_B!(workstate1::Array{Float64,1},s.H::SparseMatrixCSC{Float64,Int64},workobs1::Array{Float64,1})
	iBx[:]=iBx[:]+workstate1[:]
	#iBx[:]=iBx[:]+s.H'*(s.R \ (s.H * x))
	
    return iBx
end

# LocalWords:  iB divand

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
