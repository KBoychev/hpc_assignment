#include "functions.h"


void solve(int r,double *K,double *u,double *F){

	
	
	int dgesv_piv[r];
	int dgesv_inf;
	int dgesv_c=1;
	char dgesv_t='T';

	dgetrf_(&r, &r, K, &r, dgesv_piv, &dgesv_inf);
	dgetrs_(&dgesv_t, &r, &dgesv_c, K, &r, dgesv_piv, F, &r, &dgesv_inf);

	

}