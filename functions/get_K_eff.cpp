///////////////////////////////////////////////////////////////////////////////////
// get_K_eff() Populates the elements of the effective stiffness matrix [K_eff]. 
// ------------------------------------------------------------------------------
// @param dt <double> - Stepsize (s)
// @param M <double*> - Mass matrix 
// @param n <int> - Degrees of freedom
// @param K <double*> - Stiffness matrix (Banded symmetric storage)
// @param N_n <int> - Number of nodes
// ------------------------------------------------------------------------------

#include "functions.h"


void get_K_eff(double &dt, double *M, int &n, double *K, int &N_n) {

	for(int n_n=0;n_n<N_n;n_n++){

		K[(3*n_n)*5+4]+=4.0/(dt*dt)*M[3*n_n];
		K[(3*n_n+1)*5+4]+=4.0/(dt*dt)*M[3*n_n+1];
		K[(3*n_n+2)*5+4]+=4.0/(dt*dt)*M[3*n_n+2];

	}

}