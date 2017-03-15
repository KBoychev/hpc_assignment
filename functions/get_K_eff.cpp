#include "functions.h"


void get_K_eff(double &dt, double *M, int n_k, double *K) {

	for(int j=0;j<n_k;j++){
		K[4*n_k+j]+=4.0/(dt*dt)*M[j];
	}

}