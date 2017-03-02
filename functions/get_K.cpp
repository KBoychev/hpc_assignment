#include "functions.h"


void get_K(int r, int c, double *K, int r_e, double *K_e, int N) {

	int n = 1;

	for (int i = 0; i<r; i++) {

		for (int j = i; j<c; j++) {

			K[i*r + j] = 0;

			if (i >= 3 * n && i <= (3 * n + 2) && n<N) {

				if (j <= (3 * (n + 1) + 2)) {

					if (j <= (3 * n + 2)) {

						K[i*r + j] = K_e[((i - 3 * (n - 1)) - 3)*r_e + ((j - 3 * (n - 1)) - 3)] + K_e[(i - 3 * (n - 1))*r_e + (j - 3 * (n - 1))];

					}
					else {

						K[i*r + j] = K_e[((i - 3 * (n - 1)) - 3)*r_e + ((j - 3 * (n - 1)) - 3)];

					}

				}


			}
			else {

				if (j <= (3 * n + 2)) {
					K[i*r + j] = K_e[(i - 3 * (n - 1))*r_e + (j - 3 * (n - 1))];
				}


			}

			if (j != i) {
				K[j*r + i] = K[i*r + j];
			}

		}

		if (i == (3 * n + 2)) {
			n++;
		}

	}

	for (int i = 0; i<r; i++) {

		for (int j = 0; j<c; j++) {
			if(i==0||i==1||i==2||i==r-3||i==r-2||i==r-1){
				if(j==i){
					K[i*r+j]=1;
				}else{
					K[i*r+j]=0;				
				}
			}
		}
	}

}