#include "functions.h"


void get_F(double &t, double &T, double &l, double *F, int &N_n) {

	for(int n_n=0;n_n<N_n;n_n++){

		F[3*n_n+0]=0;

		if(n_n==(N_n-1)/2){			
			F[3*n_n+1]=t*1000.0/T*(l+1);
		}else{
			F[3*n_n+1]=t*1000.0/T*l;
		}

		F[3*n_n+2]=0;

	}

}