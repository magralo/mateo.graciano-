#include "bayesm.h"

//EXTRA FUNCTIONS SPECIFIC TO THE MAIN FUNCTION--------------------------------------------
vec drawwi(vec const& w, vec const& mu, mat const& sigmai, int p, int y){

// Wayne Taylor 9/8/2014

//function to draw w_i by Gibbing thru p vector

  int above;
	double bound;
  vec outwi = w;
  vec maxInd(2);

	for(int i = 0; i<p; i++){	
		bound = 0.0;
		for(int j = 0; j<p; j++) if(j!=i) {
        maxInd[0] = bound;
        maxInd[1] = outwi[j];
        bound = max(maxInd);}
    
    if (y==(i+1))
			above = 0;
		else 
			above = 1;
    
		vec CMout = condmom(outwi,mu,sigmai,p,i+1);
    // outwi[i] = rtrun1(CMout[0],CMout[1],bound,above);
    outwi[i] = trunNorm(CMout[0],CMout[1],bound,above);
  }

  return (outwi);
}

vec draww(vec const& w, vec const& mu, mat const& sigmai, ivec const& y){

// Wayne Taylor 9/8/2014 

//function to gibbs down entire w vector for all n obs
  
  int n = y.n_rows;
  int p = sigmai.n_cols;
  int ind; 
  vec outw = zeros<vec>(w.n_rows);
  
	for(int i = 0; i<n; i++){
    ind = p*i;
		outw.subvec(ind,ind+p-1) = drawwi(w.subvec(ind,ind+p-1),mu.subvec(ind,ind+p-1),sigmai,p,y[i]);
	}

  return (outw);
}