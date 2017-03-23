///////////////////////////////////////////////////////////////////////////////////
// Populates the elements of the stiffness matrix [K]. 
// ------------------------------------------------------------------------------
// @param E <double> - Young's modulus (Pa)
// @param A <double> - Cross-section area (m^2)
// @param I <double> - Second moment of area (m^4)
// @param l <double> - Element length (m)
// @param n <int> - Degrees of freedom
// @param K <double*> - Stiffness matrix (Banded symmetric storage)
// @param N_n <int> - Number of nodes
// ------------------------------------------------------------------------------


#include "functions.h"


void get_K(double &E, double &A,  double &I, double &l, int &n, double *K, int &N_n) {


	

	for(int n_n=0;n_n<N_n;n_n++){	

		if(n_n<N_n-1){

			K[0*n+(3*n_n+5)]=(6.0*E*I)/(l*l);

			K[1*n+(3*n_n+3)]=-(A*E)/l;
			K[1*n+(3*n_n+4)]=-(12.0*E*I)/(l*l*l);
			K[1*n+(3*n_n+5)]=(2.0*E*I) / l;
			
			K[2*n+(3*n_n+4)]=-(6.0*E*I)/(l*l);		
		}


		K[4*n+(3*n_n)]=2.0*(A*E)/l;
		K[4*n+(3*n_n+1)]=2.0*(12.0*E*I)/(l*l*l);
		K[4*n+(3*n_n+2)]=2.0*(4.0*E*I)/l;
	}


}