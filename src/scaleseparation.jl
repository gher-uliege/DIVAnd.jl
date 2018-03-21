"""

    Phi1,H1Phi1,Phi2,H2Phi2=scaleseparation(K1andH1K1,K2andH2K2,d;niter=10)    

# Input: 

* `K1andH1K1` : function, when called with a vector of data d, provides in return Kd,d-HK d ie the gridded analysis Kd an residual d-HKd for analysis tool 1

* `K2andH2K2` : function, when called with a vector of data d, provides in return Kd,d-HK d ie the gridded analysis Kd an residual d-HKd for analysis tool 2

* `d` : data array 

* `niter=` : optional keyword parameter defining the number of iterations used to invers I - H2K2 H1K2. Default 10

# Output: 
      	
* `Phi1` : analysis for tool 1
* `H1Phi1`  : analysis at data locations for tool 1
* `Phi2` : analysis for tool 2
* `H2Phi2`  : analysis at data locations for tool 2
 	
Tool to separate scales using two different analysis 	
K1 should be related to the larger scales (or scales with high signal/noise ratios)
see 	Multi-scale optimal interpolation: application to DINEOF analysis spiced with a local optimal interpolation
	http://hdl.handle.net/2268/165394
the two fields can have different supports. Only the observational operators must provide the same data array
	
	
"""


function scaleseparation(K1andH1K1,K2andH2K2,d;niter=10)



    #Allows for different state vector size; In total 3+2*niter analyses
	# K1 should be related to the larger scales (or scales with high signal/noise ratios)
	# see 	Multi-scale optimal interpolation: application to DINEOF analysis spiced with a local optimal interpolation
	# 	http://hdl.handle.net/2268/165394
	# the two fields can have different supports. Only the observational operators must provide the same data array
	
    Phi1,w1=K1andH1K1(d)
    H1Phi1=deepcopy(w1)
    w1=d-w1
    w2=deepcopy(w1)    
    # Loop
    for i=1:niter
        bidon2,w2=K2andH2K2(w2)
        bidon1,w2=K1andH1K1(w2)
        w2=w1+w2
    end
    Phi2,w1=K2andH2K2(w2)
    H2Phi2=deepcopy(w1)
    bidon1,w1=K1andH1K1(w1)
    Phi1=Phi1-bidon1
    H1Phi1=H1Phi1-w1
    return Phi1,H1Phi1,Phi2,H2Phi2
    
end