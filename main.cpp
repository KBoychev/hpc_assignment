

#include <iostream>
#include <iomanip>
#include <string>

#include "functions.h"


int main(int argc, char *argv[]) {


	int N_e = 0;
	double L = 0;	
	double A = 0;
	double I = 0;
	double E = 0;
	double rho = 0;
	double T = 0;
	int N_t = 0;
	int eq;
	int sch;
	
	//// Beam and element geometry and material properties
	//-----------------------------------------------------
	if (argc>16){

		for (int i=1; i<argc; i++){


			if(std::strcmp(argv[i],"-L")==0){
				try{
					L=std::atof(argv[i+1]);
				}catch(const std::invalid_argument&){
					std::cout << "L argument is invalid!\n";
					return 1;
				}catch(const std::out_of_range&){
					std::cout << "L argument is out of range for a double!\n";
					return 1;
				}		
			}
			if(std::strcmp(argv[i],"-N_e")==0){
				try{
					N_e=std::stoi(argv[i+1]);
				}catch(const std::invalid_argument&){
					std::cout << "N_e argument is invalid!\n";
					return 1;
				}catch(const std::out_of_range&){
					std::cout << "N_e argument is out of range for a int!\n";
					return 1;
				}
			}
			if(std::strcmp(argv[i],"-A")==0){
				try{
					A=std::stod(argv[i+1]);
				}catch(const std::invalid_argument&){
					std::cout << "A argument is invalid!\n";
					return 1;
				}catch(const std::out_of_range&){
					std::cout << "A argument is out of range for a double!\n";
					return 1;
				}	A=std::stod(argv[i+1]);
			}
			if(std::strcmp(argv[i],"-I")==0){
				try{
					I=std::stod(argv[i+1]);
				}catch(const std::invalid_argument&){
					std::cout << "I argument is invalid!\n";
					return 1;
				}catch(const std::out_of_range&){
					std::cout << "I argument is out of range for a double!\n";
					return 1;
				}	
			}
			if(std::strcmp(argv[i],"-E")==0){
				try{
					E=std::stod(argv[i+1]);
				}catch(const std::invalid_argument&){
					std::cout << "E argument is invalid!\n";
					return 1;
				}catch(const std::out_of_range&){
					std::cout << "E argument is out of range for a double!\n";
					return 1;
				}	
			}
			if(std::strcmp(argv[i],"-T")==0){
			 	try{
					T=std::stod(argv[i+1]);
				}catch(const std::invalid_argument&){
					std::cout << "T argument is invalid!\n";
					return 1;
				}catch(const std::out_of_range&){
					std::cout << "T argument is out of range for a double!\n";
					return 1;
				}
			}
			if(std::strcmp(argv[i],"-N_t")==0){
			 	try{
					N_t=std::stoi(argv[i+1]);
				}catch(const std::invalid_argument&){
					std::cout << "N_t argument is invalid!\n";
					return 1;
				}catch(const std::out_of_range&){
					std::cout << "N_t argument is out of range for a int!\n";
					return 1;
				}
			}
			if(std::strcmp(argv[i],"-rho")==0){
			 	try{
					rho=std::stod(argv[i+1]);
				}catch(const std::invalid_argument&){
					std::cout << "Rho argument is invalid!\n";
					return 1;
				}catch(const std::out_of_range&){
					std::cout << "Rho argument is out of range for a double!\n";
					return 1;
				}
			}
			if(std::strcmp(argv[i],"-eq")==0){
				try{
					eq=std::stoi(argv[i+1]);
				}catch(const std::invalid_argument&){
					std::cout << "Eq argument is invalid!\n";
					return 1;
				}catch(const std::out_of_range&){
					std::cout << "Eq argument is out of range for a int!\n";
					return 1;
				}
			}

			if(std::strcmp(argv[i],"-sch")==0){
				try{
					sch=std::stoi(argv[i+1]);
				}catch(const std::invalid_argument&){
					std::cout << "Sch argument is invalid!\n";
					return 1;
				}catch(const std::out_of_range&){
					std::cout << "Sch argument is out of range for a int!\n";
					return 1;
				}
			}
		}

	}else{
		std::cout << "Not enough arguments!\n";
		return 1;
	}

	double l = L / N_e;
	double dt = T / ( N_t - 1);

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
	
	//// Set element mass matrix [M_e]
	//--------------------------------------------------------

	set_M_e(r_e,c_e, M_e, rho, A, l);

	disp(r_e,c_e,M_e,"M_e");

	//// Set element stiffness matrix [K_e]
	//--------------------------------------------------------

	set_K_e(r_e, c_e,K_e, E, A, I, l);

	disp(r_e, c_e, K_e, "K_e");


	//// Get global mass matrix [M]
	//--------------------------------------------------------

	get_M(r, c, M, r_e, M_e, N_e);

	disp(r,c,M,"M");

	//// Get global stiffness matrix [K] and set boundary conditions
	//--------------------------------------------------------

	get_K(r, c, K, r_e, K_e, N_e);


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

	disp(r, c, K, "K");


	//// Solve for static problem [K]{u}={F}. Solution
	// 	 is done with conjugate gradient descent algorithm,
	//   implemented with BLAS library routines.
	//--------------------------------------------------------

	std::cout<<eq<<std::endl;
	std::cout<<sch<<std::endl;

	if(eq==0){



		double qy=1000;
		double Fy=1000;

		//// Set element load vector {F_e}
		//--------------------------------------------------------

		set_F_e(r_e, F_e, l,qy);

		disp(r_e,1,F_e,"F_e");

		//// Get global load vector {F} and set boundary conditions
		//--------------------------------------------------------

		get_F(r, F, r_e, F_e, N_e,Fy);

		F[(r-1)/2]=F[(r-1)/2]+Fy;

		F[0]=0;
		F[1]=0;
		F[2]=0;

		F[r-3]=0;
		F[r-2]=0;
		F[r-1]=0;

		disp(r,1,F,"F");


		std::cout<<std::endl;
	 	std::cout<<"Solving [K]{u}={F}"<<std::endl;

		inv(K,r); //The inverse function ovverwrites the matrix [K] !

		disp(r,c,K,"K_inv");

		for(int i=0;i<r;i++){
			u[i]=0;
			for(int j=0;j<c;j++){
				u[i]=u[i]+K[i*r+j]*F[j];
			}
		}

		std::cout<<"Done!"<<std::endl;

	 	std::cout<<std::endl;

		disp(r,1,u,"u");

	}

	if(eq==1){

	//// Solve the dynamic problem [M]d2{u}/dt2+[K]{u}={F}. Solution
	//	 is performed with explicit integration method.
	//--------------------------------------------------------

		std::cout<<"Solving [M]d2{u}/dt2+[K]{u}={F}"<<std::endl;
		std::cout<<std::endl;


		double t=0;
		double qy=0;
		double Fy=0;


		double *u = new double[r]();
		double *u_n = new double[r](); 
		double *tmp = new double[r]();

		if(sch==0){

			
			double *u_p = new double[r](); 
			


			for(int n_t=0;n_t<=N_t;n_t++){

				std::cout<<"Iteration "<<n_t<<" t="<<t<<std::endl;

				qy=t*1000.0/T;

				Fy=t*1000.0/T;

				set_F_e(r_e, F_e, l,qy);

				get_F(r, F, r_e, F_e, N_e,Fy);

				F[(r-1)/2]=F[(r-1)/2]+Fy;

				F[0]=0;
				F[1]=0;
				F[2]=0;

				F[r-3]=0;
				F[r-2]=0;
				F[r-1]=0;

				//Get u(n+1)
				//-------------------------------
				for(int i=0;i<r;i++){

					tmp[i]=0;

					for(int j=0;j<c;j++){

						tmp[i]=tmp[i]-(K[i*r+j]-2.0/(dt*dt)*M[i*r+j])*u[j]-1.0/(dt*dt)*M[i*r+j]*u_p[j];

					}

					tmp[i]=tmp[i]+F[i];

					u_n[i]=tmp[i]*dt*dt*1.0/M[i*r+i];
				}

				//Set boundary conditions
				//-------------------------------

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

		}


		if(sch==1){

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

			for(int n_t=0;n_t<=N_t;n_t++){

				std::cout<<"Iteration "<<n_t<<" t="<<t<<std::endl;

				qy=t*1000.0/T;

				Fy=t*1000.0/T;

				set_F_e(r_e, F_e, l,qy);

				get_F(r, F, r_e, F_e, N_e,Fy);

				F[(r-1)/2]=F[(r-1)/2]+Fy;

				F[0]=0;
				F[1]=0;
				F[2]=0;

				F[r-3]=0;
				F[r-2]=0;
				F[r-1]=0;

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

				t=t+dt;

			}

		}
	

		std::cout<<std::endl;

		std::cout<<"Done!"<<std::endl;

	 	std::cout<<std::endl;

		disp(r,1,u,"u");

	}

	return 0;
}
