#include "functions.h"


void get_M(int &r, int &c, double *M, int &r_e, double *M_e, int &N) {

	int n = 1;
	int i_e=0;
	int j_e=0;

	for(int n=1;n<=N;n++){

		i_e=0;

		for(int i=(3*n-3);i<3*(n+1);i++){

			j_e=0;

			for(int j=(3*n-3);j<3*(n+1);j++){

					M[i*r+j]=M_e[i_e*r_e+j_e]+M[i*r+j];

					j_e++;
			}

			i_e++;
		}

	}

}