

#include <iostream>
#include <iomanip>


#include "functions.h"


int main() {



	//// Beam and element geometry and material properties
	//--------------------------------------------------------

	int N_e = 4; // - (number of elements)
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

	//// Free memory from element matrices that are not going to be used again
	//--------------------------------------------------------

	delete[] M_e;
	delete[] K_e;

	//// Solve for static problem [K]{u}={F}. Solution
	// 	 is done with conjugate gradient descent algorithm,
	//   implemented with BLAS library routines.
	//--------------------------------------------------------

	

	// cout<<endl;
 	// 	cout<<"Solving [K]{u}={F}"<<endl;

	// inv(K,r); //The inverse function ovverwrites the matrix [K] !

	// disp(r,c,K,"K_inv");

	// for(int i=0;i<r;i++){
	// 	u[i]=0;
	// 	for(int j=0;j<c;j++){
	// 		u[i]=u[i]+K[i*r+j]*F[j];
	// 	}
	// }

	// cout<<"Done!"<<endl;

 	// cout<<endl;

	// disp(r,1,u,"u");



	//// Solve the dynamic problem [M]d2{u}/dt2+[K]{u}={F}. Solution
	//	 is performed with explicit integration method.
	//--------------------------------------------------------

	// cout<<"Solving [M]d2{u}/dt2+[K]{u}={F}"<<endl;

	// cout<<endl;


	 double *tmp = new double[r]();

	
	// for(int n_t=1;n_t<=N_t;n_t++){

	// 	//cout<<"Iteration "<<n_t<<" t="<<t<<endl;

	// 	qy=t*1000.0/T;

	// 	Fy=t*1000.0/T;

	// 	set_Fe(r_e, F_e, l,qy);

	// 	get_F(r, F, r_e, F_e, N_e,Fy);


	// 	//Get u(n+1)
	// 	//-------------------------------
	// 	for(int i=0;i<r;i++){

	// 		tmp[i]=0;

	// 		for(int j=0;j<c;j++){

	// 			tmp[i]=tmp[i]-(K[i*r+j]-2.0/(dt*dt)*M[i*r+j])*u[j]-1.0/(dt*dt)*M[i*r+j]*u_p[j];

	// 		}

	// 		tmp[i]=tmp[i]+F[i];

	// 		u_n[i]=tmp[i]*dt*dt*1.0/M[i*r+i];
	// 	}

	// 	//Set boundary conditions
	// 	//-------------------------------

	// 	u_n[0]=0;
	// 	u_n[1]=0;
	// 	u_n[2]=0;
	// 	u_n[r-1]=0;
	// 	u_n[r-2]=0;
	// 	u_n[r-3]=0;

		
	// 	for(int i=0;i<r;i++){
	// 		u_p[i]=u[i];
	// 		u[i]=u_n[i];
	// 	}
		

	// 	t=t+dt;


	// }


	// delete[] u_p;

	// cout<<endl;

	// cout<<"Done!"<<endl;

 	// 	cout<<endl;

	// disp(r,1,u,"u");
	
	

	// //// Solve the dynamic problem [M]d2{u}/dt2+[K]{u}={F}. Solution
	// //	 is performed with implicit integration method.
	// //--------------------------------------------------------

	cout<<"Solving [M]d2{u}/dt2+[K]{u}={F}"<<endl;

	memset(u,0,r*sizeof(*u));

	double *u_tt = new double[r]();
	double *u_tt_n = new double[r]();
	double *u_t = new double[r]();
	double *u_t_n = new double[r]();
	double *K_eff = new double[r*c]();


	for(int i=0;i<r;i++){
		for(int j=0;j<r;j++){
			K_eff[i*r+j]=4.0/(dt*dt)*M[i*r+j]+K[i*r+j];
		}
	}

	disp(r,c,K_eff,"K_eff");

	inv(K_eff,r);

	disp(r,c,K_eff,"K_eff");


	cout<<endl;
	
	for(int n_t=1;n_t<=2;n_t++){

		//cout<<"Iteration "<<n_t<<" t="<<t<<endl;

		qy=(t+dt)*1000.0/T;

		Fy=(t+dt)*1000.0/T;

		set_Fe(r_e, F_e, l,qy);

		get_F(r, F, r_e, F_e, N_e,Fy);


		//Get u_n
		//-------------------------------
		for(int i=0;i<r;i++){

			tmp[i]=0;

			for(int j=0;j<c;j++){

				tmp[i]=tmp[i]+M[i*r+j]*(4.0/(dt*dt)*u[j]+4.0/dt*u_t[j]+u_tt[j]);

			}

			tmp[i]=tmp[i]+F[i];
		}


		for(int i=0;i<r;i++){

			u_n[i]=0;

			for(int j=0;j<c;j++){
				u_n[i]=u_n[i]+K_eff[i*r+j]*tmp[j];

			}

			u_tt_n[i]=4.0/(dt*dt)*(u_n[i]-u[i])-4.0/dt*u_t[i]-u_tt[i];

			u_t_n[i]=u[i]+dt*(1-0.5)*u_tt[i]+dt*0.5*u_tt_n[i];
		}

		disp(r,1,F,"F");
		disp(r,1,tmp,"tmp");
		disp(r,1,u_n,"u");
		disp(r,1,u_tt_n,"u_tt");
		disp(r,1,u_t_n,"u_t");

		//Set boundary conditions
		//-------------------------------

		u_n[0]=0;
		u_n[1]=0;
		u_n[2]=0;
		u_n[r-1]=0;
		u_n[r-2]=0;
		u_n[r-3]=0;

		u_t_n[0]=0;
		u_t_n[1]=0;
		u_t_n[2]=0;
		u_t_n[r-1]=0;
		u_t_n[r-2]=0;
		u_t_n[r-3]=0;

		u_tt_n[0]=0;
		u_tt_n[1]=0;
		u_tt_n[2]=0;
		u_tt_n[r-1]=0;
		u_tt_n[r-2]=0;
		u_tt_n[r-3]=0;

		
		for(int i=0;i<r;i++){
			u[i]=u_n[i];
			u_t[i]=u_t_n[i];
			u_tt[i]=u_tt_n[i];
		}

		disp(r,1,u,"u");
		

		t=t+dt;

	}


	delete[] u_tt;
	delete[] u_tt_n;
	delete[] u_t;
	delete[] u_t_n;
	delete[] tmp;

	cout<<endl;

	cout<<"Done!"<<endl;

 	cout<<endl;

	disp(r,1,u,"u");


	return 0;
}
