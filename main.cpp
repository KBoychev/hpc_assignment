

#include <iostream>
#include <iomanip>
#include <string>
#include <fstream>
#include <cmath>


#include "cblas.h"
#include "mpi.h"
#include "functions.h"


using namespace std;


#define F77NAME(x) x##_
extern "C" {
	
	void F77NAME(dpbsv) (const char& uplo, const int& n, const int& kd, const int& nrhs, const double * a, const int& lda, double * b, const int& ldb, int& info);      
	void F77NAME(pdpbsv)(const char& uplo,const int& n, const int& bw,const int& nrhs, const double * a, const int& ja, const int *desca, const double *b,const int &ib,const int *descb, const double *work, const int &lwork,const int& info);        
	void F77NAME(pdgbsv)(const int& n, const int& kl,const int& ku, const int& nrhs, const double * a, const int& ja, const int *desca, const int *ipiv, const double *b,const int &ib,const int *descb, const double *work, const int &lwork,const int& info);        

	void Cblacs_get(int, int, int*);
	void Cblacs_pinfo(int*,int*);
	void Cblacs_gridinit(int*,char*,int,int);
	void Cblacs_gridinfo(int,int*,int*,int*,int*);
	void Cblacs_exit(int);
	void Cblacs_gridexit(int);  
}

int main(int argc, char *argv[]) {


	//// Set MPI variables and initialise MPI 
	//=============================================================


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
	double Tl = 0.5;
	int N_t = 0;
	double t=0;
	int eq;
	int sch;

	for (int i=1; i<argc; i++){

		if(strcmp(argv[i],"-L")==0){
			L=stod(argv[i+1]);
		}
		if(strcmp(argv[i],"-N_e")==0){
			N_e=stoi(argv[i+1]);
		}
		if(strcmp(argv[i],"-A")==0){
			A=stod(argv[i+1]);
		}
		if(strcmp(argv[i],"-I")==0){
			I=stod(argv[i+1]);	
		}
		if(strcmp(argv[i],"-E")==0){
			E=stod(argv[i+1]);	
		}
		if(strcmp(argv[i],"-T")==0){
			T=stod(argv[i+1]);
		}
		if(strcmp(argv[i],"-N_t")==0){
			N_t=stoi(argv[i+1]);
		}
		if(strcmp(argv[i],"-rho")==0){
			rho=stod(argv[i+1]);
		}
		if(strcmp(argv[i],"-eq")==0){
			eq=stoi(argv[i+1]);
		}
		if(strcmp(argv[i],"-sch")==0){
			sch=stoi(argv[i+1]);
		}
	}

	double l = L / N_e; //Element length

	double dt = T / ( N_t - 1); //Timestep


	//// Serial code
	//=============================================================

	if(MPI_N_P==1){

		int N_n=N_e-1; //Number of nodes

		int n = 3*N_n; //Dimension of global mass matrices [M], [K], and vector {F}

		//// Global mass [M], stiffness [K], and force {F} matrices 
		//	 and vector
		//--------------------------------------------------------

		double *M = new double[n](); // Global mass matrix [M] (Banded symmetric)
		double *K = new double[5*n](); // Global stiffness matrix [K] (Banded symmetric)
		double *F = new double[n](); //Global force vector {F}

		//// Get and display global mass matrix [M]
		//--------------------------------------------------------

		get_M(rho,A,l,M,N_n);

		disp(1,n,M,"M");

		//// Get and display global stiffness matrix [K]
		//--------------------------------------------------------

		get_K(E,A,I,l,n,K,N_n);

		disp(5,n, K, "K");

		
		//// Solve the static problem [K]{u}={F}
		//=============================================================

		if(eq==0){	
			
			cout<<endl;
			cout<<"Solving [K]{u}={F}"<<endl;

		 	// Global displacement vector {u}
			//------------------------------

			double *u = new double[n](); 

			// Get and display the global force vector {F} at time t=1s.
			// ------------------------------
			
			t=1;
			T=1;

			get_F(t,T,l,F,N_n);

			F[3*((N_n-1)/2)+1]=t*1000.0/T*(l+1);

			disp(n,1,F,"F");


			rmjr(5,n,K);

			// Get u 
			//------------------------------	

			int dpbsv_info=0;


			cblas_dcopy(n,F,1,u,1);

			F77NAME(dpbsv)('U',n,4,1,K,5,u,n,dpbsv_info);

			
			if (dpbsv_info) {
				cout << "Solution error "<< dpbsv_info <<" !"<< endl;
			}

			// Store and display u 
			//------------------------------

			log(N_n,l,u,"task_static.log");

			cout<<endl<<"Done!"<<endl;

			disp(n,1,u,"u");

		}

		//// Solve the dynamic problem [M]d2{u}/dt2+[K]{u}={F} on
		// 	 one process.
		//=============================================================

		if(eq==1){

				//Explicit scheme
				//=============================================================

			if(sch==0){

				

				std::ofstream center_node_log_file;

				center_node_log_file.open("center_node_defletion_wtr_time.log");

				if(center_node_log_file.good()){
					
					

					std::cout<<"Solving [M]d2{u}/dt2+[K]{u}={F} with explicit scheme"<<std::endl;
					std::cout<<std::endl;


					rmjr(5,n,K);

					double *u = new double[n]();
					double *u_p = new double[n]();
					double *u_m = new double[n]();
					double *tmp = new double[n]();

					t=0;

					for(int n_t=0;n_t<=N_t;n_t++){

						center_node_log_file<<t<<","<<u[3*((N_n-1)/2)+1]<<"\n";
						
						get_F(t,T,l,F,N_n);

						for(int n_n=0;n_n<N_n;n_n++){

							F[3*n_n+0]=0;

							if(t<=Tl){	
								F[3*n_n+1]=t*1000.0/Tl*l;	
							}else{
								F[3*n_n+1]=1000.0*l;
							}

							F[3*n_n+2]=0;
						}

						if(t<=Tl){	
							F[3*((N_n-1)/2)+1]=t*1000.0/Tl*(l+1);
						}else{
							F[3*((N_n-1)/2)+1]=1000.0*(l+1);
						}

						cblas_dsbmv(CblasColMajor,CblasUpper,n,4,-dt*dt,K,5,u,1,0.0,tmp,1);

						for(int i=0;i<n;i++){
							tmp[i]+=dt*dt*(F[i]+2.0/(dt*dt)*M[i]*u[i]-1.0/(dt*dt)*M[i]*u_m[i]);
							u_p[i]=tmp[i]/M[i];

						}

						cblas_dcopy(n,u,1,u_m,1);
						cblas_dcopy(n,u_p,1,u,1);

						t=t+dt;
					}

					center_node_log_file.close();

					log(N_n,l,u,"task_dynamic_explicit.log");

					cout<<endl<<"Done!"<<endl;

					disp(n,1,u,"u");
					
				}

			}
			

				//Implicit scheme 
				//=============================================================

			if(sch==1){


				int dpbsv_info=0;


				std::cout<<"Solving [M]d2{u}/dt2+[K]{u}={F} with implicit scheme"<<std::endl;
				std::cout<<std::endl;

				get_K_eff(dt,M,n,K);

				disp(5,n,K,"K_eff");

				rmjr(5,n,K);

				double *u = new double[n]();
				double *u_p = new double[n]();
				double *u_t = new double[n]();
				double *u_t_p = new double[n]();
				double *u_tt = new double[n]();
				double *u_tt_p = new double[n]();
				double *tmp = new double[n]();
				double *K_tmp = new double[5*n]();

				t=0;

				for(int n_t=0;n_t<=N_t;n_t++){

					t=t+dt;

					get_F(t,T,l,F,N_n);


					
					F[3*((N_n-1)/2)+1]=t*1000.0/T*(l+1);

					for(int i=0;i<n;i++){
						tmp[i]=F[i]+M[i]*(1.0/(0.25*dt*dt)*u[i]+1.0/(0.25*dt)*u_t[i]+(1.0/(2.0*0.25)-1.0)*u_tt[i]);
					}

					cblas_dcopy(n,tmp,1,u_p,1);

					cblas_dcopy(5*n,K,1,K_tmp,1);

					F77NAME(dpbsv)('U',n,4,1,K_tmp,5,u_p,n,dpbsv_info);

					if (dpbsv_info) {
						cout << "Solution error "<< dpbsv_info <<" !"<< endl;
					}

					for(int i=0;i<n;i++){
						u_tt_p[i]=1.0/(0.25*dt*dt)*(u_p[i]-u[i])-1.0/(0.25*dt)*u_t[i]-(1.0/(2.0*0.25)-1.0)*u_tt[i];  
						u_t_p[i]=u_t[i]+dt*0.5*u_tt[i]+dt*0.5*u_tt_p[i];
					}

					cblas_dcopy(n,u_p,1,u,1);
					cblas_dcopy(n,u_t_p,1,u_t,1);
					cblas_dcopy(n,u_tt_p,1,u_tt,1);

					
				}


				log(N_n,l,u,"task_dynamic_implicit.log");

				cout<<endl<<"Done!"<<endl;

				disp(n,1,u,"u");


			}

		}

	}

	//// Parallel code
	//=============================================================
	if(MPI_N_P==2){


		//// Solve the dynamic problem [M]d2{u}/dt2+[K]{u}={F} on
		// 	 two processes.
		//=============================================================

		if(eq==1){


				//Explicit scheme
				//=============================================================

			if(sch==0){

					int N_n=(N_e-2)/2+1; //Number of nodes

					int n = 3*N_n; //Dimension of global mass matrices [M], [K], and vector {F}

					//// Global mass [M], stiffness [K], and force {F} matrices 
					//	 and vector
					//--------------------------------------------------------

					double *M = new double[n](); // Global mass matrix [M] (Banded symmetric)
					double *K = new double[5*n](); // Global stiffness matrix [K] (Banded symmetric)
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
						disp(5,n,K,"K");
					}

					if(MPI_P_ID==0){
						K[4*n+(3*(N_n-1)+0)]=(A*E)/l;
						K[4*n+(3*(N_n-1)+1)]=(12.0*E*I)/(l*l*l);
						K[4*n+(3*(N_n-1)+2)]=(4.0*E*I)/l;
					}
					
					if(MPI_P_ID==1){
						K[4*n+0]=(A*E)/l;
						K[4*n+1]=(12.0*E*I)/(l*l*l);
						K[4*n+2]=(4.0*E*I)/l;
					}


					if(MPI_P_ID==0){
						std::cout<<"Solving [M]d2{u}/dt2+[K]{u}={F} with explicit scheme on two processes"<<std::endl;
						std::cout<<std::endl;
					}


					rmjr(5,n,K);

					double *MPI_BUFF = new double[3]();
					double *u = new double[n]();
					double *u_p = new double[n]();
					double *u_m = new double[n]();
					double *tmp = new double[n]();

					t=0;

					for(int n_t=0;n_t<=N_t;n_t++){

						get_F(t,T,l,F,N_n);

						if(MPI_P_ID==0){
							F[3*(N_n-1)+1]=t*1000.0/T*(l+1);
						}

						if(MPI_P_ID==1){
							F[3*0+1]=t*1000.0/T*(l+1);
						}

						cblas_dsbmv(CblasColMajor,CblasUpper,n,4,-dt*dt,K,5,u,1,0.0,tmp,1);

						if(MPI_N_P==2 && MPI_P_ID==0){

							MPI_BUFF[0]=tmp[n-3];
							MPI_BUFF[1]=tmp[n-2];
							MPI_BUFF[2]=tmp[n-1];

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

							tmp[n-3]+=MPI_BUFF[0];
							tmp[n-2]+=MPI_BUFF[1];
							tmp[n-1]+=MPI_BUFF[2];
						}

						for(int i=0;i<n;i++){
							tmp[i]+=dt*dt*(F[i]+2.0/(dt*dt)*M[i]*u[i]-1.0/(dt*dt)*M[i]*u_m[i]);
							u_p[i]=tmp[i]/M[i];
						}

						cblas_dcopy(n,u,1,u_m,1);

						cblas_dcopy(n,u_p,1,u,1);

						t=t+dt;
					}



					if(MPI_P_ID==0){

						log(N_n,l,u,"task_dynamic_explicit_parallel.log");

						cout<<endl<<"Done!"<<endl;

						disp(n,1,u,"u");

					}	
					

				}

				//Implicit scheme 
				//=============================================================

				if(sch==1){

					int N_n=(N_e-2)/2+1;

					if(MPI_P_ID==1){
						N_n--;// For scalapack odd number of nodes no overlap
					}

					int n = 3*N_n; //Dimension of global mass matrices [M], [K], and vector {F}
				
					//// Global mass [M], stiffness [K], and force {F} matrices 
					//	 and vector
					//--------------------------------------------------------

					double *M = new double[n](); // Global mass matrix [M] (Banded)
					double *K = new double[17*n](); // Global stiffness matrix [K] (Banded)
					double *F = new double[n](); //Global force vector {F}

					//// Get and display global mass matrix [M]
					//--------------------------------------------------------

					get_M(rho,A,l,M,N_n);

					if(MPI_P_ID==0){
						disp(1,n,M,"M");

					}

					//// Get and display global stiffness matrix [K]
					//--------------------------------------------------------

					

					if(MPI_P_ID==0){

						for(int n_n=0;n_n<N_n;n_n++){	
							
							if(n_n<N_n-1){

								K[(0+8)*n+(3*n_n+5)]=(6.0*E*I)/(l*l);

								K[(1+8)*n+(3*n_n+3)]=-(A*E)/l;
								K[(1+8)*n+(3*n_n+4)]=-(12.0*E*I)/(l*l*l);
								K[(1+8)*n+(3*n_n+5)]=(2.0*E*I) / l;
								
								K[(2+8)*n+(3*n_n+4)]=-(6.0*E*I)/(l*l);		
							}

							K[(6+8)*n+(3*n_n+2)]=-(6.0*E*I)/(l*l);

							K[(7+8)*n+(3*n_n+0)]=-(A*E)/l;
							K[(7+8)*n+(3*n_n+1)]=-(12.0*E*I)/(l*l*l);
							K[(7+8)*n+(3*n_n+2)]=(2.0*E*I) / l;
							
							K[(8+8)*n+(3*n_n+1)]=(6.0*E*I)/(l*l);		


							K[(4+8)*n+(3*n_n)]=2*(A*E)/l;
							K[(4+8)*n+(3*n_n+1)]=2*(12.0*E*I)/(l*l*l);
							K[(4+8)*n+(3*n_n+2)]=2*(4.0*E*I)/l;
						}
					}

					if(MPI_P_ID==1){

						for(int n_n=0;n_n<N_n;n_n++){	
							

							K[(0+8)*n+(3*n_n+2)]=(6.0*E*I)/(l*l);

							K[(1+8)*n+(3*n_n+0)]=-(A*E)/l;
							K[(1+8)*n+(3*n_n+1)]=-(12.0*E*I)/(l*l*l);
							K[(1+8)*n+(3*n_n+2)]=(2.0*E*I) / l;
							
							K[(2+8)*n+(3*n_n+1)]=-(6.0*E*I)/(l*l);		
							
							if(n_n<N_n-1){

								K[(6+8)*n+(3*n_n+2)]=-(6.0*E*I)/(l*l);

								K[(7+8)*n+(3*n_n+0)]=-(A*E)/l;
								K[(7+8)*n+(3*n_n+1)]=-(12.0*E*I)/(l*l*l);
								K[(7+8)*n+(3*n_n+2)]=(2.0*E*I) / l;
								
								K[(8+8)*n+(3*n_n+1)]=(6.0*E*I)/(l*l);


							}

							K[(4+8)*n+(3*n_n)]=2*(A*E)/l;
							K[(4+8)*n+(3*n_n+1)]=2*(12.0*E*I)/(l*l*l);
							K[(4+8)*n+(3*n_n+2)]=2*(4.0*E*I)/l;
						}
					}




					if(MPI_P_ID==0){
						disp(17,n,K,"K");
					}

					

					//Scalapack setup 
					//-------------------------------

					char order='R';
					int mype, npe, ctx, myrow, mycol;
					int nrow=1;
					int ncol=2;
					
					Cblacs_pinfo(&mype,&npe);
					Cblacs_get(0,0,&ctx);
					Cblacs_gridinit(&ctx,&order,1,npe);
					Cblacs_gridinfo(ctx,&nrow,&ncol,&myrow,&mycol);

					int desc_K[7];
					desc_K[0]=501; //Type is banded matrix 1-by-P (LHS)
					desc_K[1]=ctx; //Context
					desc_K[2]=3*(N_e-1); //Problem size
					desc_K[3]=3*((N_e-2)/2+1); //Blocking
					desc_K[4]=0; //Process row/column
					desc_K[5]=17; //Local size
					desc_K[6]=0; //Reserved

					int desc_u_p[7]; 
					desc_u_p[0]=502; //Type is banded matrix P-by-1 (RHS)
					desc_u_p[1]=ctx; //Context
					desc_u_p[2]=3*(N_e-1); //Problem size
					desc_u_p[3]=3*((N_e-2)/2+1); //Blocking
					desc_u_p[4]=0; //Process row/column
					desc_u_p[5]=3*((N_e-2)/2+1); //Local size
					desc_u_p[6]=0; //Reserved

					
					int *ipiv=new int[3*(N_e-1)]; //Ipiv
					int lwork=(3*((N_e-2)/2+1)+4)*(4+4)+6*(4+4)*(4+2*4)+max(1*(3*((N_e-2)/2+1)+2*4+4*4),1); // Work size
					double *work=new double[lwork]; // Work

					int info; //Info

					if(MPI_P_ID==0){
						std::cout<<"Solving [M]d2{u}/dt2+[K]{u}={F} with implicit scheme on two processes"<<std::endl;
						std::cout<<std::endl;
					}

					for(int j=0;j<n;j++){
						K[(4+8)*n+j]+=4.0/(dt*dt)*M[j];
					}

					if(MPI_P_ID==0){
						disp(17,n,K,"K_eff");
					}

					rmjr(17,n,K);

					double *u = new double[n]();
					double *u_p = new double[n]();
					double *u_t = new double[n]();
					double *u_t_p = new double[n]();
					double *u_tt = new double[n]();
					double *u_tt_p = new double[n]();
					double *tmp = new double[n]();
					double *K_tmp = new double[17*n]();

					t=0;

					for(int n_t=0;n_t<=N_t;n_t++){

						t=t+dt;

						get_F(t,T,l,F,N_n);

						if(MPI_P_ID==0){
							F[3*(N_n-1)+1]=t*1000.0/T*(l+1);
						}

						
						for(int i=0;i<n;i++){
							tmp[i]=F[i]+M[i]*(1.0/(0.25*dt*dt)*u[i]+1.0/(0.25*dt)*u_t[i]+(1.0/(2.0*0.25)-1.0)*u_tt[i]);
						}

						cblas_dcopy(n,tmp,1,u_p,1);

						cblas_dcopy(17*n,K,1,K_tmp,1);


						//Solve the system K_tml * u_p'=u_p (i.e. RHS vector replaced by soltion)

						//Problem size
						//Number of subdiagonals
						//Number of superdiagonals
						//Number of RHS
						//Matrix
						//Matrix index
						//Matrix descriptor
						//Ipiv
						//Vector
						//Vector index
						//Vector descriptor
						//Work
						//Size of work
						//Info
						
						F77NAME(pdgbsv)(3*(N_e-1),4,4,1,K_tmp,1,desc_K,ipiv,u_p,1,desc_u_p,work,lwork,info);

						if (info) {
							cout << "Solution error "<< info <<" !"<< endl;
						}  
						

						for(int i=0;i<n;i++){
							u_tt_p[i]=1.0/(0.25*dt*dt)*(u_p[i]-u[i])-1.0/(0.25*dt)*u_t[i]-(1.0/(2.0*0.25)-1.0)*u_tt[i];  
							u_t_p[i]=u_t[i]+dt*0.5*u_tt[i]+dt*0.5*u_tt_p[i];
						}

						cblas_dcopy(n,u_p,1,u,1);
						cblas_dcopy(n,u_t_p,1,u_t,1);
						cblas_dcopy(n,u_tt_p,1,u_tt,1);
						
					}

					Cblacs_gridexit(ctx);

					if(MPI_P_ID==0){

						log(N_n,l,u,"task_dynamic_implicit_parallel.log");

						cout<<endl<<"Done!"<<endl;

						disp(n,1,u,"u");

					}

				}



			}

		}


		MPI_Finalize();

		return 0;
	}
