///////////////////////////////////////////////////////////////////////////////////
// get_F() Populates the elements of the force vector {F}. 
// ------------------------------------------------------------------------------
// @param t <double> - Time (s)
// @param T <double> - Period (s)
// @param l <double> - Element length (m)
// @param F <double*> - Force vector 
// @param N_n <int> - Number of nodes
// ------------------------------------------------------------------------------

#include "functions.h"


void get_F(double &t, double &Tl, double &l, double *F, int &N_n) {

	for(int n_n=0;n_n<N_n;n_n++){

		F[3*n_n+0]=0;

		if(t<Tl){
			F[3*n_n+1]=t*1000.0/Tl*l;
		}else{
			F[3*n_n+1]=1000.0*l;
		}

		F[3*n_n+2]=0;

	}

}
