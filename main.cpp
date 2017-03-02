

#include <iostream>
#include <iomanip>

#include <Accelerate/Accelerate.h>


#include "functions.h"


int main() {



	//// Beam and element geometry and material properties
	//--------------------------------------------------------

	int N_e = 24; // - (number of elements)
	double L = 10.0; // m (beam length)
	double l = L / N_e; //m (element length)

	double A = 0.1*0.12; //m^2 (cross-sectional area)
	double I = (0.1*0.12*0.12*0.12)/12; //m^4 (second moment of area)
	double E = 210000000000; //2.1*pow(10,11); //Pa (Young's modulus)
	double rho = 7850; //kg/m^3 (material density)
	double T = 1; //s (end time)
	int N_t = 10000; //- (number of time steps)
	double dt=T/(N_t-1); //s (timestep)
	double t=0; //s (time)
	double qy=1000;
	double Fy=1000;


	cout<<setprecision(6);

	cout<<"l:"<<l<<"m"<<endl;
	cout<<"A:"<<A<<"m2"<<endl;
	cout<<"I:"<<I<<"m4"<<endl;
	cout<<"E:"<<E<<"Pa"<<endl;
	cout<<"dt:"<<dt<<"s"<<endl;


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

	get_M(r, c, M, r_e, M_e, N_e);

	disp(r,c,M,"M");

	//// Get global stiffness matrix [K] and set boundary conditions
	//--------------------------------------------------------

	get_K(r, c, K, r_e, K_e, N_e);

	disp(r, c, K, "K");

	//// Get global load vector {F} and set boundary conditions
	//--------------------------------------------------------

	get_F(r, F, r_e, F_e, N_e,Fy);

	disp(r,1,F,"F");


	//// Solve for static problem [K]{u}={F}. Solution
	// 	 is done with conjugate gradient descent algorithm,
	//   implemented with BLAS library routines.
	//--------------------------------------------------------

	//fucks up [K] matrix for the next problem!!!!! separate codes

	// cout<<endl;
 	// cout<<"Solving [K]{u}={F}"<<endl;

	// int dgesv_piv[r];
	// int dgesv_inf;
	// int dgesv_c=1;
	// char dgesv_t='T';

	// dgetrf_(&r, &r, K, &r, dgesv_piv, &dgesv_inf);
	// dgetrs_(&dgesv_t, &r, &dgesv_c, K, &r, dgesv_piv, F, &r, &dgesv_inf);

	// cout<<"Done!"<<endl;

 	// 	cout<<endl;

	// disp(r,1,F,"u");


	cout<<"Solving [M]d2{u}/dt2+[K]{u}={F}"<<endl;

	cout<<endl;

	//// Solve the dynamic problem. [M]d2{u}/dt2+[K]{u}={F}
	//--------------------------------------------------------
	
	for(int n_t=1;n_t<=N_t;n_t++){

		//cout<<"Iteration "<<n_t<<" t="<<t<<endl;

		qy=t*1000.0/T;

		Fy=t*1000.0/T;

		set_Fe(r_e, F_e, l,qy);

		get_F(r, F, r_e, F_e, N_e,Fy);


		//Get u(n+1)
		//-------------------------------
		for(int i=0;i<r;i++){

			u_n[i]=0;

			for(int j=0;j<c;j++){

				u_n[i]=u_n[i]-(K[i*r+j]-2.0/(dt*dt)*M[i*r+j])*u[j]-1.0/(dt*dt)*M[i*r+j]*u_p[j];

			}

			u_n[i]=u_n[i]+F[i];
			u_n[i]=u_n[i]*dt*dt*1.0/M[i*r+i];
		}

		//Set boundary conditions

		u_n[0]=0;
		u_n[1]=0;
		u_n[2]=0;
		u_n[r-1]=0;
		u_n[r-2]=0;
		u_n[r-3]=0;

		
		for(int i=0;i<r;i++){
			u_p[i]=u[i];
			u[i]=u_n[i];
		}
		

		t=t+dt;


	}

	cout<<endl;

	cout<<"Done!"<<endl;

 	cout<<endl;

	disp(r,1,u,"u");
	
	

	return 0;
}
