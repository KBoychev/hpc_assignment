#include "functions.h"


void get_K(int &r, int &c, double *K, int &r_e, double *K_e, int &N) {


	int i_e=0;
	int j_e=0;

	for(int n=1;n<=N;n++){

		i_e=0;

		for(int i=(3*n-3);i<3*(n+1);i++){

			j_e=0;

			for(int j=(3*n-3);j<3*(n+1);j++){

					K[i*r+j]=K_e[i_e*r_e+j_e]+K[i*r+j];

					j_e++;
			}

			i_e++;
		}

	}

}