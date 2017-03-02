#include "functions.h"


void set_Fe(int n_e, double *F_e, double l, double qy) {

	F_e[0] = 0;
	F_e[1] = qy / 2;
	F_e[2] = (qy * l) / 12;
	F_e[3] = 0;
	F_e[4] = qy / 2;
	F_e[5] = -(qy * l) / 12;

}
