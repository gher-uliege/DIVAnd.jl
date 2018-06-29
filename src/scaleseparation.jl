"""

    Phi1,H1Phi1,Phi2,H2Phi2=scaleseparation(K1andH1K1,K2andH2K2,d;niter=10)    

# Input: 

* `K1andH1K1` : function, when called with a vector of data d, provides in return Kd,HK d ; i.e. the gridded analysis Kd and analysis HKd at data points for analysis tool 1

* `K2andH2K2` : function, when called with a vector of data d, provides in return Kd,HK d ; i.e. the gridded analysis Kd and analysis HKd at data points for analysis tool 2

* `d` : data array 

* `niter=` : optional keyword parameter defining the number of iterations used to invert I - H2K2 H1K2. Default is 10

# Output: 
      	
* `Phi1` : analysis for tool 1 in which analysis of scale 2 is taken out
* `H1Phi1`  : analysis at data locations for tool 1 in which analysis of scale 2 is taken out
* `Phi2` : analysis for tool 2 in which analysis of scale 1 is taken out
* `H2Phi2`  : analysis at data locations for tool 2 in which analysis of scale 1 is taken out
 	
Tool to separate scales using two different analysis provided as two input functions

K1 should be related to the larger scales (or scales with high signal/noise ratios) and K2 to smaller or less energetic scales. If in doubt invert both and test
with different number of iterations while looking at convergence.

see 	"Multi-scale optimal interpolation: application to DINEOF analysis spiced with a local optimal interpolation"
	http://hdl.handle.net/2268/165394
	
	
Here the two fields can have different supports (one could be a 3D analysis and the other one a season-depth analysis for example. Only the observational operators must provide the same data array at the output. In other words K1,HK1=K1andH1K1 should provide an output array HK1 of the same dimensions as the data array d and the output HK2 from K2,HK2=K2andH2K2
	
	
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
