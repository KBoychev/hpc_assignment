#include "functions.h"


void get_M_p(int p,int &r, int &c, double *M, double *M_p) {

	
	int c_p=(r+1)/2+1;
	int r_p=(r+1)/2+1;

	for(int i=0;i<r_p;i++){

		for(int j=0;j<c_p;j++){

			M_p[i*c_p+j]=M[(i+p*6)*c+(j+p*6)];
				
		}

	}

}