#include "functions.h"


void get_F_p(int p,int &r, double *F, int &r_p, double *F_p) {

	for(int i=0;i<r_p;i++){
			F_p[i]=F[(i+p*6)];
	}
}