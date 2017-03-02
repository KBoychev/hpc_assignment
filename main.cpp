

#include <iostream>
#include <iomanip>
#include <string>
#include <cmath>
#include <Accelerate/Accelerate.h>


using namespace std;



void disp(int r,int c, double *m , string l) {

	cout << endl;

	cout << setprecision(6);

	if (c == 1) {
		cout << r << "x" << c << " {" << l << "}=" << endl;
		for (int i = 0; i < r; i++) {
			cout << setw(12);
			cout << m[i];
			cout << endl;
		}
	}
	else {
		cout << r << "x" << c << " [" << l << "]=" << endl;

		for (int i = 0; i < r; i++) {
			for (int j = 0; j < c; j++) {
				cout << setw(12);
				cout << m[i*r + j];
			}
			cout << endl;
		}
	}

	

	

	cout << endl;

}


void solve(int r,double *K,double *u,double *F){

	
	
	int dgesv_piv[r];
	int dgesv_inf;
	int dgesv_c=1;
	char dgesv_t='T';

	dgetrf_(&r, &r, K, &r, dgesv_piv, &dgesv_inf);
	dgetrs_(&dgesv_t, &r, &dgesv_c, K, &r, dgesv_piv, F, &r, &dgesv_inf);

	

}

void set_Me(int r_e, int c_e, double *Me, double rho, double A, double l) {

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

 
void set_K_e(int r_e, int c_e, double *K_e, double E, double A, double I, double l) {

	for (int i = 0; i<r_e; i++) {
		for (int j = i; j<c_e; j++) {

			K_e[i*r_e + j] = 0;

			if (j == i) {

				if(i==0||i==3){
					K_e[i*r_e + j] = (A*E) / l;
				}

				if (i == 2 || i == 5) {

					K_e[i*r_e + j] = (4.0*E*I) / l;
				}

				if (i == 1 || i == 4) {
					K_e[i*r_e + j] = (12.0*E*I) / (l*l*l);
				}

			}

			if (i == 0 && j == 3) {
				K_e[i*r_e + j] = -(A*E) / l;
			}
			if (i == 1 && j == 2) {
				K_e[i*r_e + j] = (6.0*E*I) / (l*l);
			}

			if (i == 4 && j == 5) {
				K_e[i*r_e + j] = -(6.0*E*I) / (l*l);
			}
			if (i == 1 && j == 4) {
				K_e[i*r_e + j] = -(12.0*E*I) / (l*l*l);
			}
			if (i == 1 && j == 5) {
				K_e[i*r_e + j] = (6.0*E*I) / (l*l);
			}

			if (i == 2 && j == 4) {
				K_e[i*r_e + j] = -(6.0*E*I) / (l*l);
			}
			if (i == 2 && j == 5) {
				K_e[i*r_e + j] = (2.0*E*I) / l;

			}

			if (j != i) {
				K_e[j*r_e + i] = K_e[i*r_e + j];
			}
		}
	}
}

void set_Fe(int n_e, double *F_e, double l, double qy) {

	F_e[0] = 0;
	F_e[1] = qy / 2;
	F_e[2] = (qy * l) / 12;
	F_e[3] = 0;
	F_e[4] = qy / 2;
	F_e[5] = -(qy * l) / 12;

}

void get_M(int r, int c, double *M, int r_e, double *M_e, int N_e) {

	int n_e = 1;

	for (int i = 0; i<r; i++) {

		if (i >= 3 * n_e && i <= (3 * n_e + 2) && n_e<N_e) {

			M[i*r+i]=M_e[((i - 3 * (n_e - 1)) - 3)*r_e + ((i - 3 * (n_e - 1)) - 3)] + M_e[(i - 3 * (n_e - 1))*r_e + (i - 3 * (n_e - 1))];

		}else{

			M[i*r + i] = M_e[(i - 3 * (n_e - 1))*r_e + (i - 3 * (n_e - 1))];
				
		}

	}


}

void set_bc_K(int r,int c,double *K){

	for (int i = 0; i<r; i++) {

		for (int j = 0; j<c; j++) {

				if(i==0||i==1||i==2||i==r-3||i==r-2||i==r-1){
					if(j==i){
						K[i*r+j]=1;
					}else{
						K[i*r+j]=0;
						
					}

				}


		}

	}

}

void get_K(int r, int c, double *K, int r_e, double *K_e, int N) {

	int n = 1;

	for (int i = 0; i<r; i++) {

		for (int j = i; j<c; j++) {

			K[i*r + j] = 0;

			if (i >= 3 * n && i <= (3 * n + 2) && n<N) {

				if (j <= (3 * (n + 1) + 2)) {

					if (j <= (3 * n + 2)) {

						K[i*r + j] = K_e[((i - 3 * (n - 1)) - 3)*r_e + ((j - 3 * (n - 1)) - 3)] + K_e[(i - 3 * (n - 1))*r_e + (j - 3 * (n - 1))];

					}
					else {

						K[i*r + j] = K_e[((i - 3 * (n - 1)) - 3)*r_e + ((j - 3 * (n - 1)) - 3)];

					}

				}


			}
			else {

				if (j <= (3 * n + 2)) {
					K[i*r + j] = K_e[(i - 3 * (n - 1))*r_e + (j - 3 * (n - 1))];
				}


			}

			if (j != i) {
				K[j*r + i] = K[i*r + j];
			}

		}

		if (i == (3 * n + 2)) {
			n++;
		}

	}

}


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

int main() {



	//// Beam and element geometry and material properties
	//--------------------------------------------------------

	int N_e = 8; // - (number of elements)
	double L = 10.0; // m (beam length)
	double l = L / N_e; //m (element length)

	double A = 0.1*0.12; //m^2 (cross-sectional area)
	double I = (0.1*pow(0.12,3))/12; //m^4 (second moment of area)
	double E = 2.1*pow(10,11); //Pa (Young's modulus)
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

	//// Get global stiffness matrix [M]
	//--------------------------------------------------------

	get_K(r, c, K, r_e, K_e, N_e);

	disp(r, c, K, "K");

	//// Get global load vector {F}
	//--------------------------------------------------------

	get_F(r, F, r_e, F_e, N_e,Fy);

	disp(r,1,F,"F");

	//// Set boundary conditions on [K]. At first and last node
	// no displacements in x or y.
	//-------------------------------------------------------- 

	cout<<"Setting boundary conditions on [K]"<<endl;

	set_bc_K(r,c,K);

	disp(r, c, K, "K");

	//// Solve for static equilibrium [K]{u}={F}. Solution
	// 	 is done with conjugate gradient descent algorithm,
	//   implemented with BLAS library routines.
	//--------------------------------------------------------

	cout<<endl;
 	cout<<"Solving [K]{u}={F}"<<endl;

	solve(r,K,u,F);

	cout<<"Done!"<<endl;

 	cout<<endl;

	disp(r,1,F,"u");


	//// Solve the dynamic problem.
	//--------------------------------------------------------


	// cblas_dscal(r,0,u,1); //reset u


	// cout<<endl;
 	 // 	cout<<"Solving [M]{u_tt}+[K]{u}={F} for N_t "<<N_t<<" and dt "<<dt<<endl;

	// for (int n_t=0;n_t<=N_t;n_t++){

	// 	cout<<endl;
	// 	cout<<"Iteration "<<n_t<<" t "<<t<<endl;

		   
	// 	qy=t*1000/T;
	// 	Fy=t*1000/T;

	// 	set_Fe(r_e, F_e, l,qy);

	// 	get_F(r, F, r_e, F_e, N_e,Fy);

	// 	for(int i=0;i<r;i++){

	// 		for(int j=0;j<c;j++){

	// 			u_n[i]=u_n[i]-(K[i*r+j]-2.0/(dt*dt)*M[i*r+j])*u[j]-1.0/(dt*dt)*M[i*r+j]*u_p[j];

	// 		}

	// 		u_n[i]=u_n[i]+F[i];

	// 		//cout<<u_n[i]<<endl;

	// 		u_n[i]=u_n[i]*dt*dt/M[i*r+i];

	// 		cout<<u_n[i]<<endl;

			

	// 	}

	// 	disp(r,1,u_n,"u_n");

	// 	cblas_dcopy(r,u,1,u_p,1);
	// 	cblas_dcopy(r,u_n,1,u,1);

		

	// 	t=t+dt;
	// }

	return 0;
}
