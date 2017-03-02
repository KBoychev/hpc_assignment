#include "functions.h"


void get_M(int r, int c, double *M, int r_e, double *M_e, int N_e) {

	int n_e = 1;

	for (int i = 0; i<r; i++) {

		if (i >= 3 * n_e && i <= (3 * n_e + 2) && n_e<N_e) {

			M[i*r+i]=M_e[((i - 3 * (n_e - 1)) - 3)*r_e + ((i - 3 * (n_e - 1)) - 3)] + M_e[(i - 3 * (n_e - 1))*r_e + (i - 3 * (n_e - 1))];

		}else{

			M[i*r + i] = M_e[(i - 3 * (n_e - 1))*r_e + (i - 3 * (n_e - 1))];
				
		}

	}


}