

#include <iostream>
#include <iomanip>
#include <string>


#include <Accelerate/Accelerate.h>


#include "mpi.h"
#include "functions.h"


int main(int argc, char *argv[]) {


	//// Set MPI variables and initialise MPI (TODO)
	//-----------------------------------------------------

	int MPI_N_P;
	int MPI_P_ID;
	int MPI_status;

	MPI_status = MPI_Init(&argc, &argv);

	MPI_status = MPI_Comm_size(MPI_COMM_WORLD,&MPI_N_P);

 	MPI_status = MPI_Comm_rank(MPI_COMM_WORLD,&MPI_P_ID);

	int N_e = 0;
	double L = 0;	
	double A = 0;
	double I = 0;
	double E = 0;
	double rho = 0;
	double T = 0;
	int N_t = 0;
	double t=0;
	int eq;
	int sch;

	for (int i=1; i<argc; i++){

		if(std::strcmp(argv[i],"-L")==0){
			L=std::atof(argv[i+1]);
		}
		if(std::strcmp(argv[i],"-N_e")==0){
			N_e=std::stoi(argv[i+1]);
		}
		if(std::strcmp(argv[i],"-A")==0){
			A=std::stod(argv[i+1]);
		}
		if(std::strcmp(argv[i],"-I")==0){
			I=std::stod(argv[i+1]);	
		}
		if(std::strcmp(argv[i],"-E")==0){
			E=std::stod(argv[i+1]);	
		}
		if(std::strcmp(argv[i],"-T")==0){
			T=std::stod(argv[i+1]);
		}
		if(std::strcmp(argv[i],"-N_t")==0){
			N_t=std::stoi(argv[i+1]);
		}
		if(std::strcmp(argv[i],"-rho")==0){
			rho=std::stod(argv[i+1]);
		}
		if(std::strcmp(argv[i],"-eq")==0){
			eq=std::stoi(argv[i+1]);
		}
		if(std::strcmp(argv[i],"-sch")==0){
			sch=std::stoi(argv[i+1]);
		}
	}

	double l = L / N_e; //element length
	double dt = T / ( N_t - 1); //timestep
	
	int N_n=N_e-1; //Number of nodes

	int n = 3*N_n; //Dimension of global mass matrices [M], [K], and vector {F}

	//// Global mass [M], stiffness [K], and force {F} matrices 
	//	 and vector
	//--------------------------------------------------------

	double *M = new double[n](); // Global mass matrix [M] (Banded)
	double *K = new double[9*n](); // Global stiffness matrix [K] (Banded)
	double *F = new double[n](); //Global force vector {F}

	//// Get and display global mass matrix [M]
	//--------------------------------------------------------

	get_M(rho,A,l,M,N_n);

	if(MPI_P_ID==0){

		disp(1,n,M,"M");

	}

	//// Get and display global stiffness matrix [K]
	//--------------------------------------------------------

	get_K(E,A,I,l,n,K,N_n);

	if(MPI_P_ID==0){
		disp(9,n, K, "K");
	}

	//// Solve the static problem [K]{u}={F}
	//--------------------------------------------------------

	if(eq==0){	
			
		std::cout<<std::endl;
	 	std::cout<<"Solving [K]{u}={F}"<<std::endl;

	 	// Global displacement vector {u}
		//------------------------------

	 	double *u = new double[n](); 

		// Get and display the global force vector {F} at time t=1s.
		// ------------------------------
		
		t=1;
		T=1;

		get_F(t,T,l,F,N_n);

		disp(n,1,F,"F");

		// Get u
		//------------------------------	

		char dpbsv_uplo='U';
		int dpbsv_kd=4;
		int dpbsv_nrhs=1;
		int dpbsv_ldk=5;
		int dpbsv_info;

		dpbsv_(&dpbsv_uplo,&n,&dpbsv_kd,&dpbsv_nrhs,K,&dpbsv_ldk,F,&n,&dpbsv_info);

		std::cout<<dpbsv_info<<std::endl;

		// Store and display u 
		//------------------------------

		log(n,u,"task_static.log");

 		std::cout<<std::endl<<"Done!"<<std::endl;

		disp(n,1,u,"u");

	}

	//// Solve the dynamic problem [M]d2{u}/dt2+[K]{u}={F} on
	// 	 one or two processes.
	//--------------------------------------------------------

	if(eq==1){

		int n_p;

		if(MPI_N_P==1){
			n_p=3*N_n;
		}

		if(MPI_N_P==2){
			n_p=3*((N_n-1)/2+1);
		}

		std::cout<<n_p<<std::endl;


		double *MPI_BUFF = new double[3](); //buffer of {u} to keep values of displacement for one node

		double *u; // {u} at timelayer n
		double *u_p; // {u} at timelayer n+1
		double *u_p_k;
		double *u_m; // {u} at timelayer n+1
		double *u_tt; // second time derivative of {u} at timelayer n
		double *u_tt_p; // second time derivative of {u} at timelayer n+1
		double *u_t; // first time derivative of {u} at timelayer n
		double *u_t_p; // first time derivative of {u} at timelayer n+1
		double *tmp = new double[n_p]();
		double *b;
		double R=0;
		int k=0;


		if(sch==0){

			delete u_p_k;
			delete u_tt;
			delete u_tt_p;
			delete u_t;
			delete u_t_p;
			delete b;

			u=new double[n_p]();
			u_p=new double[n_p]();
			u_m=new double[n_p]();


		}

		if(sch==1){

			delete u_m;

			u=new double[n_p]();
			u_p=new double[n_p]();
			u_p_k=new double[n_p]();
			u_tt=new double[n_p]();
			u_tt_p=new double[n_p]();
			u_t=new double[n_p]();
			u_t_p=new double[n_p]();
			b=new double[n_p]();

			get_K_eff(dt,M,n,K);	

			if(MPI_P_ID==0){
				disp(9,n,K,"K_eff");
			}

		}


		if(MPI_P_ID==0){
			std::cout<<"Solving [M]d2{u}/dt2+[K]{u}={F}"<<std::endl;
			std::cout<<std::endl;
		}
		

		//Loop for N_t timesteps
		//-------------------------------

		t=0;

		for(int n_t=0;n_t<=N_t;n_t++){


			//Explicit scheme
			//-------------------------------

			if(sch==0){

				get_F(t,T,l,F,N_n);
				
				for(int j=0;j<n_p;j++){

					tmp[j]=0;

					for(int i=0;i<9;i++){

						if(j<4){
							if(i>(3-j) && i!=4){
								tmp[j]-=K[i*n+j]*u[i-(4-j)];
							}
						}

						if(j>=4 && j<(n_p-4)){
							if(i!=4){
								tmp[j]-=K[i*n+j]*u[i+(j-4)];
							}

						}

						if(j>=(n_p-4)){
							if(i<(4+(n_p-j))&& i!=4){
								tmp[j]-=K[i*n+j]*u[n_p-(4+(n_p-j))+i];
							}
						}
					}			
				}

				

				//mpi stuff;

				if(MPI_N_P==2 && MPI_P_ID==0){

					MPI_BUFF[0]=tmp[n_p-3];
					MPI_BUFF[1]=tmp[n_p-2];
					MPI_BUFF[2]=tmp[n_p-1];

					MPI_Send(MPI_BUFF,3,MPI_DOUBLE,1,0,MPI_COMM_WORLD);
				}

				if(MPI_N_P==2 && MPI_P_ID==1){

					MPI_BUFF[0]=tmp[0];
					MPI_BUFF[1]=tmp[1];
					MPI_BUFF[2]=tmp[2];

					MPI_Send(MPI_BUFF,3,MPI_DOUBLE,0,0,MPI_COMM_WORLD);				
				}

				if(MPI_N_P==2 && MPI_P_ID==1){

					MPI_Recv(MPI_BUFF,3,MPI_DOUBLE,0,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);

					tmp[0]+=MPI_BUFF[0];
					tmp[1]+=MPI_BUFF[1];
					tmp[2]+=MPI_BUFF[2];	

				}
			
				if(MPI_N_P==2 && MPI_P_ID==0){

					MPI_Recv(MPI_BUFF,3,MPI_DOUBLE,1,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);

					tmp[n_p-3]+=MPI_BUFF[0];
					tmp[n_p-2]+=MPI_BUFF[1];
					tmp[n_p-1]+=MPI_BUFF[2];
				}


				for(int i=0;i<n_p;i++){

					tmp[i]+=-K[4*n+i]*u[i]; //diagonal 

					if(MPI_N_P==2 && MPI_P_ID==0){
						tmp[i]+=F[i]+2.0/(dt*dt)*M[i]*u[i]-1.0/(dt*dt)*M[i]*u_m[i]; 
					}

					if(MPI_N_P==2 && MPI_P_ID==1){
						tmp[i]+=F[3*(N_e/2-1)+i]+2.0/(dt*dt)*M[i]*u[i]-1.0/(dt*dt)*M[i]*u_m[i];
					}

				    if(MPI_N_P==1){
				    	
				    	tmp[i]+=F[i]+2.0/(dt*dt)*M[i]*u[i]-1.0/(dt*dt)*M[i]*u_m[i];
				    }
			
					u_p[i]=tmp[i]*dt*dt*1.0/M[i];
				}


				for(int i=0;i<n_p;i++){
					u_m[i]=u[i];					
					u[i]=u_p[i];
				}

				t=t+dt;
			}


			//Implicit scheme (Jacobi iterations)
			//-------------------------------
			if(sch==1){

				t=t+dt;

				get_F(t,T,l,F,N_n);

				k=0;
				R=1;

				while(R>1*std::pow(10,-12)){

					for(int j=0;j<n_p;j++){

						tmp[j]=0;

						for(int i=0;i<9;i++){

							if(j<4){
								if(i>(3-j) && i!=4){
									tmp[j]+=K[i*n+j]*u_p[i-(4-j)];
								}
							}

							if(j>=4 && j<(n-4)){
								if(i!=4){
									tmp[j]+=K[i*n+j]*u_p[i+(j-4)];
								}
							}

							if(j>=(n-4)){
								if(i<(4+(n-j)) && i!=4){
									tmp[j]+=K[i*n+j]*u_p[n_p-(4+(n_p-j))+i];
								}
							}
						}
						b[j]=F[j]+M[j]*(1.0/(0.25*dt*dt)*u[j]+1.0/(0.25*dt)*u_t[j]+(1.0/(2.0*0.25)-1.0)*u_tt[j]);
					}


					for(int i=0;i<n_p;i++){
						u_p_k[i]=1/K[4*n+i]*(b[i]-tmp[i]);
					}
					
					R=0;

					for(int i=0;i<n_p;i++){
						R+=u_p_k[i]-u_p[i];					
					}

					R=std::abs(R);

					for(int i=0;i<n_p;i++){
						u_p[i]=u_p_k[i];
					}

					k++;
				}

				for(int i=0;i<n_p;i++){
					u_tt_p[i]=1.0/(0.25*dt*dt)*(u_p[i]-u[i])-1.0/(0.25*dt)*u_t[i]-(1.0/(2.0*0.25)-1.0)*u_tt[i];  
    				u_t_p[i]=u_t[i]+dt*0.5*u_tt[i]+dt*0.5*u_tt_p[i];
				}


				for(int i=0;i<n_p;i++){
					u[i]=u_p[i];
					u_t[i]=u_t_p[i];
					u_tt[i]=u_tt_p[i];
				}

			}


		}


		//Display u and write u to file 
		//-------------------------------
		if(MPI_P_ID==0){

			log(n,u,"task_dynamic.log");
			
			std::cout<<std::endl<<"Done!"<<std::endl;

			disp(n_p,1,u,"u");
			
		}

	}


	MPI_Finalize();

	return 0;
}
