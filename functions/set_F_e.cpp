#include "functions.h"


void set_Fe(int n_e, double *F_e, double l, double qy) {

	F_e[0] = 0;
	F_e[1] = qy*l / 2;
	F_e[2] = (qy * l*l) / 12;
	F_e[3] = 0;
	F_e[4] = qy*l / 2;
	F_e[5] = -(qy * l*l) / 12;

}
