#include "functions.h"


void set_K_e(int r_e, int c_e, double *K_e, double E, double A, double I, double l) {

	for (int i = 0; i<r_e; i++) {
		for (int j = i; j<c_e; j++) {

			K_e[i*r_e + j] = 0;

			if (j == i) {

				if(i==0||i==3){
					K_e[i*r_e + j] = (A*E) / l;
				}

				if (i == 2 || i == 5) {

					K_e[i*r_e + j] = (4.0*E*I) / l;
				}

				if (i == 1 || i == 4) {
					K_e[i*r_e + j] = (12.0*E*I) / (l*l*l);
				}

			}

			if (i == 0 && j == 3) {
				K_e[i*r_e + j] = -(A*E) / l;
			}
			if (i == 1 && j == 2) {
				K_e[i*r_e + j] = (6.0*E*I) / (l*l);
			}

			if (i == 4 && j == 5) {
				K_e[i*r_e + j] = -(6.0*E*I) / (l*l);
			}
			if (i == 1 && j == 4) {
				K_e[i*r_e + j] = -(12.0*E*I) / (l*l*l);
			}
			if (i == 1 && j == 5) {
				K_e[i*r_e + j] = (6.0*E*I) / (l*l);
			}

			if (i == 2 && j == 4) {
				K_e[i*r_e + j] = -(6.0*E*I) / (l*l);
			}
			if (i == 2 && j == 5) {
				K_e[i*r_e + j] = (2.0*E*I) / l;

			}

			if (j != i) {
				K_e[j*r_e + i] = K_e[i*r_e + j];
			}
		}
	}
}
