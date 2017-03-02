#include "functions.h"


void get_F(int r, double *F, int r_e, double *F_e, int N,double Fy) {

	int n = 1;

	for (int i = 0; i<r; i++) {

		F[i] = 0;

		if (i >= 3 * n && i <= (3 * n + 2) && n<N) {
			F[i] = F_e[(i - 3 * (n - 1))] + F_e[(i - 3 * (n - 1)) - 3];
		}
		else {
			F[i] = F_e[(i - 3 * (n - 1))];
		}

		if (i == (3 * n + 2)) {
			n++;
		}

	}

	F[(r-1)/2]=F[(r-1)/2]+Fy;
	F[0]=0;
	F[1]=0;
	F[2]=0;
	F[r-3]=0;
	F[r-2]=0;
	F[r-1]=0;
}