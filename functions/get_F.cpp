#include "functions.h"


void get_F(int &r, double *F, int &r_e, double *F_e, int &N,double &Fy) {


	for(int i=0;i<r;i++){
		F[i]=0;
	}

	int n = 1;
	int i_e=0;

	for(int n=1;n<=N;n++){

		i_e=0;

		for(int i=(3*n-3);i<3*(n+1);i++){

			F[i]=F_e[i_e]+F[i];

			i_e++;
		}

	}

}