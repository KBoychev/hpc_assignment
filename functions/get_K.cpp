
#include "functions.h"


void get_K(double &E, double &A,  double &I, double &l, int &n, double *K, int &N_n) {


	

	for(int n_n=0;n_n<N_n;n_n++){	

		if(n_n<N_n-1){

			K[0*n+(3*n_n+5)]=(6.0*E*I)/(l*l);

			K[1*n+(3*n_n+3)]=-(A*E)/l;
			K[1*n+(3*n_n+4)]=-(12.0*E*I)/(l*l*l);
			K[1*n+(3*n_n+5)]=(2.0*E*I) / l;
			
			K[2*n+(3*n_n+4)]=-(6.0*E*I)/(l*l);		
		}


		K[4*n+(3*n_n)]=2.0*(A*E)/l;
		K[4*n+(3*n_n+1)]=2.0*(12.0*E*I)/(l*l*l);
		K[4*n+(3*n_n+2)]=2.0*(4.0*E*I)/l;
	}


}