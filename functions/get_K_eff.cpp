#include "functions.h"


void get_K_eff(double &dt, double *M, int &n, double *K) {

	for(int j=0;j<n;j++){
		K[4*n+j]+=4.0/(dt*dt)*M[j];
	}

}