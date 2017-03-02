

#include <iostream>
#include <iomanip>

#include <Accelerate/Accelerate.h>


#include "functions.h"


int main() {



	//// Beam and element geometry and material properties
	//--------------------------------------------------------

	int N_e = 2; // - (number of elements)
	double L = 10.0; // m (beam length)
	double l = L / N_e; //m (element length)

	double A = 0.1*0.12; //m^2 (cross-sectional area)
	double I = (0.1*0.12*0.12*0.12)/12; //m^4 (second moment of area)
	double E = 210000000000; //2.1*pow(10,11); //Pa (Young's modulus)
	double rho = 7850; //kg/m^3 (material density)
	double T = 1; //s (end time)
	int N_t = 10; //- (number of time steps)
	double dt=T/N_t; //s (timestep)
	double t=0; //s (time)
	double qy=1000;
	double Fy=1000;


	//// Element mass,stiffness, force matrices 
	//	 and vectors
	//--------------------------------------------------------

	int r_e = 6; //rows 
	int c_e = 6; //columns

	double *M_e = new double[r_e*c_e]();
	double *K_e = new double[r_e*c_e]();
	double *F_e = new double[r_e]();


	//// Global mass, stiffness, load, displacement matrices 
	//	 and vectors
	//--------------------------------------------------------
	int r = 3 * (N_e + 1); //rows 
	int c = 3 * (N_e + 1); //columns

	double *M = new double[r*c]();
	double *K = new double[r*c]();
	double *F = new double[r]();
	double *u = new double[r]();
	double *u_p = new double[r](); //u at previous time step
	double *u_n = new double[r](); //u at next timestep


	//// Set element mass matrix [M_e]
	//--------------------------------------------------------

	set_Me(r_e,c_e, M_e, rho, A, l);

	disp(r_e,c_e,M_e,"M_e");

	//// Set element stiffness matrix [K_e]
	//--------------------------------------------------------

	set_K_e(r_e, c_e,K_e, E, A, I, l);

	disp(r_e, c_e, K_e, "K_e");

	//// Set element load vector {F_e}
	//--------------------------------------------------------

	set_Fe(r_e, F_e, l,qy);

	disp(r_e,1,F_e,"F_e");

	//// Get global mass matrix [M]
	//--------------------------------------------------------

	get_K(r, c, M, r_e, M_e, N_e);

	disp(r,c,M,"M");

	//// Get global stiffness matrix [K] and set boundary conditions
	//--------------------------------------------------------

	get_K(r, c, K, r_e, K_e, N_e);

	disp(r, c, K, "K");

	//// Get global load vector {F} and set boundary conditions
	//--------------------------------------------------------

	get_F(r, F, r_e, F_e, N_e,Fy);

	disp(r,1,F,"F");


	//// Solve for static equilibrium [K]{u}={F}. Solution
	// 	 is done with conjugate gradient descent algorithm,
	//   implemented with BLAS library routines.
	//--------------------------------------------------------

	cout<<endl;
 	cout<<"Solving [K]{u}={F}"<<endl;

	int dgesv_piv[r];
	int dgesv_inf;
	int dgesv_c=1;
	char dgesv_t='T';

	dgetrf_(&r, &r, K, &r, dgesv_piv, &dgesv_inf);
	dgetrs_(&dgesv_t, &r, &dgesv_c, K, &r, dgesv_piv, F, &r, &dgesv_inf);

	cout<<"Done!"<<endl;

 	cout<<endl;

	disp(r,1,F,"u");


	//// Solve the dynamic problem.
	//--------------------------------------------------------


	

	return 0;
}
