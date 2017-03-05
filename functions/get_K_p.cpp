#include "functions.h"


void get_K_p(int p,int &r, int &c, double *K, double *K_p) {

	int c_p=(r+1)/2+1;
	int r_p=(r+1)/2+1;

	for(int i=0;i<r_p;i++){

		for(int j=0;j<c_p;j++){

			K_p[i*c_p+j]=K[(i+p*6)*c+(j+p*6)];
				
		}

	}

}