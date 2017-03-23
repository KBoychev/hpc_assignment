///////////////////////////////////////////////////////////////////////////////////
// get_K_eff() Populates the elements of the effective stiffness matrix [K_eff]. 
// ------------------------------------------------------------------------------
// @param dt <double> - Stepsize (s)
// @param M <double*> - Mass matrix 
// @param n <int> - Degrees of freedom
// @param K <double*> - Stiffness matrix (Banded symmetric storage)
// ------------------------------------------------------------------------------

#include "functions.h"


void get_K_eff(double &dt, double *M, int &n, double *K) {

	for(int j=0;j<n;j++){
		K[4*n+j]+=4.0/(dt*dt)*M[j];
	}

}