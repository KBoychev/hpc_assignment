#include "functions.h"


void set_Me(int &r_e, int &c_e, double *Me, double &rho, double &A, double &l) {

	for (int i = 0; i<r_e; i++) {
		for (int j = i; j<c_e; j++) {

			Me[i*r_e + j] = 0;

			if (i == j) {
				if (i == 2 || i == 5) {
					Me[i*r_e + j] = rho*A*l*1.0 / 24.0*l*l;
				}
				else {
					Me[i*r_e + j] = rho*A*l*1.0 / 2.0;
				}

			}

		}
	}
}