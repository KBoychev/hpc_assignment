#include "functions.h"


void get_M(double &rho, double &A, double &l, double *M, int &N_n) {

	for(int n_n=0;n_n<N_n;n_n++){
		
		M[3*n_n]=rho*A*l*1.0;
		M[3*n_n+1]=rho*A*l*1.0;
		M[3*n_n+2]=rho*A*l*l*l*1.0/12.0;
	}

}