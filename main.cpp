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

	void F77NAME(dpbsv)(const char& uplo, const int& n, const int& kd, const int& nrhs, const double * a, const int& lda, double * b, const int& ldb, int& info);      
	
	void F77NAME(dpbtrf)(const char& uplo, const int& n, const int& kd, const double * a, const int& lda, int& info);     

	void F77NAME(dpbtrs)(const char& uplo, const int& n, const int& kd, const int& nrhs, const double * a, const int& lda, double * b, const int& ldb, int& info); 

	void F77NAME(pdpbtrf)(const char& uplo,const int& n, const int& bw,const double * a, const int& ja, const int *desca, const double *af,const int &laf,const double *work, const int &lwork,const int& info);        

	void F77NAME(pdpbtrs)(const char& uplo,const int& n, const int& bw,const int& nrhs, const double * a, const int& ja, const int *desca, const double *b,const int &ib,const int *descb, const double *af,const int &laf,const double *work, const int &lwork,const int& info);        

	void Cblacs_get(int, int, int*);

	void Cblacs_pinfo(int*,int*);

	void Cblacs_gridinit(int*,char*,int,int);

	void Cblacs_gridinfo(int,int*,int*,int*,int*);

	void Cblacs_exit(int);

	void Cblacs_gridexit(int);  
}

int main(int argc, char *argv[]) {

	///////////////////////////////////////////////////////////////////////////////////
	// MPI setup
	//---------------------------------------------------------------------------------

	int MPI_N_P; //Number of processes
	int MPI_P_ID; // Process ID
	int MPI_status; //Status

	MPI_status = MPI_Init(&argc, &argv);

	MPI_status = MPI_Comm_size(MPI_COMM_WORLD,&MPI_N_P);

	MPI_status = MPI_Comm_rank(MPI_COMM_WORLD,&MPI_P_ID);

	///////////////////////////////////////////////////////////////////////////////////
	// Command line arguments
	//---------------------------------------------------------------------------------

	int N_e = 0; //Number of elements
	double L = 0; //Beam length
	double A = 0; //Area
	double I = 0; //Second moment of area
	double E = 0; //Young's modulus
	double rho = 0; //Density
	double T = 0; //Period
	double Tl = 0; //Loading period
	int N_t = 0; //Number of iterations
	double t=0; //Time
	int eq=0; //Equation (0-static, 1-dynamic)
	int sch=0; //Scheme (0-explicit, 1-implicit)

	for (int i=1; i<argc; i++){

		if(strcmp(argv[i],"-L")==0){
			L=stod(argv[i+1]);
		}
		if(strcmp(argv[i],"-N_e")==0){
			N_e=stoi(argv[i+1]);
			if(N_e%2!=0){
				std::cout<<"Please use an even number of elements!"<<std::endl;
				MPI_Finalize();
				return 0;
			}
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
		if(strcmp(argv[i],"-Tl")==0){
			Tl=stod(argv[i+1]);
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

	///////////////////////////////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////
	// Serial code
	//---------------------------------------------------------------------------------

	if(MPI_N_P==1){

		int N_n=N_e-1; //Number of nodes

		int n = 3*N_n; //Dimension of global mass matrices [M], [K], and vector {F}

		// IT IS NOT EFFICIENT TO CREATE THE ELEMENT MATRICES FIRST AND THEN THE GLOBAL MATRICES
		// GENERATE GLOBAL MATRICES WITHOUT USING ELEMENT MATRICES

		// Mass [M], stiffness [K], and force {F} matrices and vector

		double *M = new double[n](); // Mass matrix [M] (Banded symmetric)
		double *K = new double[5*n](); // Stiffness matrix [K] (Banded symmetric)
		double *F = new double[n](); //Force vector {F}

		// Get and display stiffness matrix [K]

		get_K(E,A,I,l,n,K,N_n);

		disp(5,n, K, "K");

		///////////////////////////////////////////////////////////////////////////////////
		///////////////////////////////////////////////////////////////////////////////////
		// Solve the static problem [K]{u}={F}
		//---------------------------------------------------------------------------------

		if(eq==0){	
			
			cout<<endl;
			cout<<"Solving [K]{u}={F}"<<endl;

		 	//Variables needed for solution
			double *u = new double[n](); //u
			int info=0;

			//Get F at t=1 and T=1 (this will give qy=1000 and Fy=1000)

			t=1;
			T=1;

			get_F(t,T,l,F,N_n);

			//Set center node loading
			F[3*((N_n-1)/2)+1]=t*1000.0/T*(l+1);

			//Display F
			disp(n,1,F,"F");

			//Convert K from row major to column major (for BLAS and LAPACK libraries)
			rm2cm(5,n,K);
						
			//Copy F into u
			cblas_dcopy(n,F,1,u,1);

			//Solve the system [K]*{u}'={u}
			F77NAME(dpbsv)('U',n,4,1,K,5,u,n,info);

			if (info) {
				cout << "DPBSV error "<< info <<" !"<< endl;
			}

			//Save results to file and display them
			log(N_n,l,u,"results/task_static.log");

			cout<<endl<<"Done!"<<endl;

			disp(n,1,u,"u");

		}

		///////////////////////////////////////////////////////////////////////////////////
		///////////////////////////////////////////////////////////////////////////////////
		// Solve the dynamic problem [M]d2{u}/dt2+[K]{u}={F} on one process.
		//---------------------------------------------------------------------------------

		if(eq==1){

			// Get and display mass matrix [M]

			get_M(rho,A,l,M,N_n);

			disp(1,n,M,"M");

			///////////////////////////////////////////////////////////////////////////////////
			// Explicit scheme
			//---------------------------------------------------------------------------------

			if(sch==0){

				//Open a new file to store the center node y displacement with respect to time
				std::ofstream center_node_log_file;

				center_node_log_file.open("results/explicit_center_node_defletion_wtr_time.log");

				if(center_node_log_file.good()){
					
					
					//Convert K from row major to column major (for BLAS and LAPACK libraries)
					rm2cm(5,n,K);

					//Variables needed for solution
					double *u = new double[n](); //u at time layer n
					double *u_p = new double[n](); //u at time layer n+1
					double *u_m = new double[n](); //u at time layer n-1
					double *tmp = new double[n](); //temporary storage used during computations

					//Initialise time to 0
					t=0;

					//Display that we're solving explicit problem
					std::cout<<"Solving [M]d2{u}/dt2+[K]{u}={F} with explicit scheme"<<std::endl;
					std::cout<<std::endl;

					//Loop for N_t iterations
					for(int n_t=0;n_t<=N_t;n_t++){

						//Save time and center node y displacement to file
						center_node_log_file<<t<<","<<u[3*((N_n-1)/2)+1]<<"\n";
						
						//Get F at iteration n
						get_F(t,Tl,l,F,N_n);

						//If t is less than Tl then change the center node loading linearly, if not keep it constant
						if(t<Tl){	
							F[3*((N_n-1)/2)+1]=t*1000.0/Tl*(l+1);
						}else{
							F[3*((N_n-1)/2)+1]=1000.0*(l+1);
						}

						//Multiply -dt*dt*[K]*u and store the result in tmp
						cblas_dsbmv(CblasColMajor,CblasUpper,n,4,-dt*dt,K,5,u,1,0.0,tmp,1);

						//Add to tmp the remaining part of the right hand side of the equation
						for(int i=0;i<n;i++){
							tmp[i]+=dt*dt*(F[i]+2.0/(dt*dt)*M[i]*u[i]-1.0/(dt*dt)*M[i]*u_m[i]);
							u_p[i]=tmp[i]/M[i]; // [M] is diagonal so the inverse is 1/[M]

						}

						//Copy u into u_m and u_p into u

						cblas_dcopy(n,u,1,u_m,1);
						cblas_dcopy(n,u_p,1,u,1);


						//Incement time
						t=t+dt;
					}

					//Close center node log file
					center_node_log_file.close();

					//Save results to file and display them
					log(N_n,l,u,"results/task_dynamic_explicit.log");

					cout<<endl<<"Done!"<<endl;

					disp(n,1,u,"u");
					
				}

			}

			///////////////////////////////////////////////////////////////////////////////////
			// Implicit scheme 
			//---------------------------------------------------------------------------------

			if(sch==1){

				//Open a new file to store the center node y displacement with respect to time
				std::ofstream center_node_log_file;

				center_node_log_file.open("results/implicit_center_node_defletion_wtr_time.log");

				if(center_node_log_file.good()){

					
					// Get and display the effective matrix [K_eff] ([K] is being overwritten with [K_eff])

					get_K_eff(dt,M,n,K);

					disp(5,n,K,"K_eff");

					//Convert K from row major to column major (for BLAS and LAPACK libraries)

					rm2cm(5,n,K);

					//Variables needed for solution
					double *u = new double[n](); //u at time layer n
					double *u_p = new double[n](); //u at time layer n+1
					double *u_t = new double[n](); //first derivative of u at time layer n
					double *u_t_p = new double[n](); //first derivative of u at time layer n+1
					double *u_tt = new double[n](); //second derivative of u at time layer n
					double *u_tt_p = new double[n](); //second derivative of u at time layer n+1
					double *tmp = new double[n](); //temporary storage used during computations
					int info=0;

					//Initialise time to 0
					t=0;

					//Display that we're solving implicit problem 
					std::cout<<"Solving [M]d2{u}/dt2+[K]{u}={F} with implicit scheme"<<std::endl;
					std::cout<<std::endl;


					F77NAME(dpbtrf)('U',n, 4, K,5,info);

					if (info) {
						cout << "DPBTRF error "<< info <<" !"<< endl;
					}

					//Loop for N_t iterations
					for(int n_t=0;n_t<=N_t;n_t++){

						//Save time and center node y displacement to file
						center_node_log_file<<t<<","<<u[3*((N_n-1)/2)+1]<<"\n";

						//Increment time
						t=t+dt;

						//Get F at iteration n+1
						get_F(t,Tl,l,F,N_n);

						//If t is less than Tl then change the center node loading linearly, if not keep it constant
						if(t<Tl){	
							F[3*((N_n-1)/2)+1]=t*1000.0/Tl*(l+1);
						}else{
							F[3*((N_n-1)/2)+1]=1000.0*(l+1);
						}

						//Calculate right hand side of equation
						for(int i=0;i<n;i++){
							tmp[i]=F[i]+M[i]*(1.0/(0.25*dt*dt)*u[i]+1.0/(0.25*dt)*u_t[i]+(1.0/(2.0*0.25)-1.0)*u_tt[i]);
						}

						//Copy the right hand side into u_p
						cblas_dcopy(n,tmp,1,u_p,1);

						//Solve the system [K]*{u_p}'={u_p}

						F77NAME(dpbtrs)('U',n,4,1,K, 5,u_p, n,info); 

						if (info) {
							cout << "DPBTRS error "<< info <<" !"<< endl;
						}

						//Calculate u_tt_p and u_t_p
						for(int i=0;i<n;i++){
							u_tt_p[i]=1.0/(0.25*dt*dt)*(u_p[i]-u[i])-1.0/(0.25*dt)*u_t[i]-(1.0/(2.0*0.25)-1.0)*u_tt[i];  
							u_t_p[i]=u_t[i]+dt*0.5*u_tt[i]+dt*0.5*u_tt_p[i];
						}

						//Copy u_p into u, u_t_p int u_t and u_tt_p int u_tt
						cblas_dcopy(n,u_p,1,u,1);
						cblas_dcopy(n,u_t_p,1,u_t,1);
						cblas_dcopy(n,u_tt_p,1,u_tt,1);

						
					}

					//Close center node log file
					center_node_log_file.close();

					//Save results to file and display them
					log(N_n,l,u,"results/task_dynamic_implicit.log");

					cout<<endl<<"Done!"<<endl;

					disp(n,1,u,"u");
				}

			}

		}

	}

	///////////////////////////////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////
	// Parallel code
	//---------------------------------------------------------------------------------

	if(MPI_N_P==2){

		///////////////////////////////////////////////////////////////////////////////////
		///////////////////////////////////////////////////////////////////////////////////
		// Solve the dynamic problem [M]d2{u}/dt2+[K]{u}={F} on two processes.
		//---------------------------------------------------------------------------------

		if(eq==1){

			///////////////////////////////////////////////////////////////////////////////////
			// Explicit scheme
			//---------------------------------------------------------------------------------

			if(sch==0){

					double *MPI_u = new double[3*(N_e-1)]; //MPI array to store the displacements from the two processes after finishing

					int N_n=(N_e-2)/2+1; //Number of nodes

					int n = 3*N_n; //Dimension of mass matrices [M], [K], and vector {F}


					// IT IS NOT EFFICIENT TO CREATE THE ELEMENT MATRICES FIRST AND THEN THE GLOBAL MATRICES
					// GENERATE GLOBAL MATRICES WITHOUT USING ELEMENT MATRICES

					// Mass [M], stiffness [K], and force {F} matrices and vector

					double *M = new double[n](); // Mass matrix [M] (Banded symmetric)
					double *K = new double[5*n](); // Stiffness matrix [K] (Banded symmetric)
					double *F = new double[n](); //Force vector {F}

					// Get and display the mass matrix [M] on process with ID 0

					get_M(rho,A,l,M,N_n);

					if(MPI_P_ID==0){

						disp(1,n,M,"M");

					}

					// Get and display the stiffness matrix [K] on process with ID 0

					get_K(E,A,I,l,n,K,N_n);

					if(MPI_P_ID==0){
						disp(5,n,K,"K");
					}

					//If process ID is 0 change the middle node properties (so we don't add them twice during multiplication exchange between processes)
					if(MPI_P_ID==0){
						K[4*n+(3*(N_n-1)+0)]=(A*E)/l;
						K[4*n+(3*(N_n-1)+1)]=(12.0*E*I)/(l*l*l);
						K[4*n+(3*(N_n-1)+2)]=(4.0*E*I)/l;
					}
					//If process ID is 1 change the middle node properties (so we don't add them twice during multiplication exchange between processes)
					if(MPI_P_ID==1){
						K[4*n+0]=(A*E)/l;
						K[4*n+1]=(12.0*E*I)/(l*l*l);
						K[4*n+2]=(4.0*E*I)/l;
					}

					
					//Convert K from row major to column major (for blas and lapack libraries)
					rm2cm(5,n,K);


					//Variables needed for solution
					double *MPI_BUFF = new double[3](); //MPI buffer for exchanging u values at the shared node (3 degrees of freedom)
					double *u = new double[n](); //u at current time layer n
					double *u_p = new double[n](); //u at time layer n+1
					double *u_m = new double[n](); //u at time layer n-1
					double *tmp = new double[n](); //temporary storage used during computations

					//Initialise time to 0
					t=0;


					//Display that we're solving explicit problem on two processes
					if(MPI_P_ID==0){
						std::cout<<"Solving [M]d2{u}/dt2+[K]{u}={F} with explicit scheme on two processes"<<std::endl;
						std::cout<<std::endl;
					}

					//Loop for N_t iterations
					for(int n_t=0;n_t<=N_t;n_t++){


						//Get the force vector {F} at time t
						get_F(t,T,l,F,N_n);


						//If process ID is 0 set last node y-force to qy*l+Fy
						if(MPI_P_ID==0){
							F[3*(N_n-1)+1]=t*1000.0/T*(l+1);
						}

						//If process ID is 1 set first node y-force to qy*l+Fy
						if(MPI_P_ID==1){
							F[3*0+1]=t*1000.0/T*(l+1);
						}

						//Multiply -dt*dt*[K]*u and store the result in tmp
						cblas_dsbmv(CblasColMajor,CblasUpper,n,4,-dt*dt,K,5,u,1,0.0,tmp,1);



						//If process ID is 0 send the multiplication results at the shared node (3 degrees of freedom) to process ID 1
						if(MPI_N_P==2 && MPI_P_ID==0){

							MPI_BUFF[0]=tmp[n-3];
							MPI_BUFF[1]=tmp[n-2];
							MPI_BUFF[2]=tmp[n-1];

							MPI_Send(MPI_BUFF,3,MPI_DOUBLE,1,0,MPI_COMM_WORLD);
						}

						//If process ID is 1 send the multiplication results at the shared node (3 degrees of freedom) to process ID 0
						if(MPI_N_P==2 && MPI_P_ID==1){

							MPI_BUFF[0]=tmp[0];
							MPI_BUFF[1]=tmp[1];
							MPI_BUFF[2]=tmp[2];

							MPI_Send(MPI_BUFF,3,MPI_DOUBLE,0,0,MPI_COMM_WORLD);				
						}

						//If process ID is 1 receive the multiplication results at the shared node (3 degrees of freedom) from process ID 0
						if(MPI_N_P==2 && MPI_P_ID==1){

							MPI_Recv(MPI_BUFF,3,MPI_DOUBLE,0,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);

							tmp[0]+=MPI_BUFF[0];
							tmp[1]+=MPI_BUFF[1];
							tmp[2]+=MPI_BUFF[2];	

						}
						
						//If process ID is 0 receive the multiplication results at the shared node (3 degrees of freedom) from process ID 1
						if(MPI_N_P==2 && MPI_P_ID==0){

							MPI_Recv(MPI_BUFF,3,MPI_DOUBLE,1,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);

							tmp[n-3]+=MPI_BUFF[0];
							tmp[n-2]+=MPI_BUFF[1];
							tmp[n-1]+=MPI_BUFF[2];
						}

						//Add to tmp the remaining part of the right hand side of the equation
						for(int i=0;i<n;i++){
							tmp[i]+=dt*dt*(F[i]+2.0/(dt*dt)*M[i]*u[i]-1.0/(dt*dt)*M[i]*u_m[i]);
							u_p[i]=tmp[i]/M[i]; // [M] is diagonal so the inverse is 1/[M]
						}

						//Copy u int u_m and u_p into u
						cblas_dcopy(n,u,1,u_m,1);
						cblas_dcopy(n,u_p,1,u,1);


						//Increment time
						t=t+dt;
					}



					//Gather u from the processes and remove overlap
					
					if(MPI_P_ID==0){

						MPI_Gather(u,n,MPI_DOUBLE,MPI_u,n,MPI_DOUBLE,0,MPI_COMM_WORLD);

					}else{

						//Remove overlapping results (i.e. remove center node results from process id 1 and gather)
						double *u_tmp=new double[n-3]();

						for(int i=0;i<n;i++){
							if(i>2){
								u_tmp[i-3]=u[i];
							}
						}
					

						MPI_Gather(u_tmp,n-3,MPI_DOUBLE,MPI_u,n-3,MPI_DOUBLE,0,MPI_COMM_WORLD);
					}


					//If process ID is 0 save results to file and display them
					if(MPI_P_ID==0){

						N_n=N_e-1;
						n=3*N_n;

						log(N_n,l,MPI_u,"results/task_dynamic_explicit_parallel.log");

						cout<<endl<<"Done!"<<endl;

						disp(n,1,MPI_u,"u");

					}	
					
				}

				///////////////////////////////////////////////////////////////////////////////////
				// Implicit scheme 
				//---------------------------------------------------------------------------------

				if(sch==1){

					double *MPI_u = new double[3*(N_e-1)];

					int N_n=(N_e-2)/2+1;

					if(MPI_P_ID==1){
						N_n--;//Process with ID 1 gets one node less
					}

					int n = 3*N_n; //Dimension of mass matrices [M], [K], and vector {F}

					// Mass [M], stiffness [K], and force {F} matrices and vector

					double *M = new double[n](); // Mass matrix [M] (Banded symmetric)
					double *K = new double[5*n](); // Stiffness matrix [K] (Banded symmetric)
					double *F = new double[n](); // Force vector {F}

					// Get and display mass matrix [M] on process with ID 0

					get_M(rho,A,l,M,N_n);

					if(MPI_P_ID==0){
						disp(1,n,M,"M");
					}

					// Get and display stiffness matrix [K] on process with ID 0

					get_K(E,A,I,l,n,K,N_n);

					
					//If process ID is 1 modify the K matrix and remove the upper left triangle of zeros and replace them with values 

					if(MPI_P_ID==1){

						K[0*n+2]=(6.0*E*I)/(l*l);

						K[1*n+0]=-(A*E)/l;
						K[1*n+1]=-(12.0*E*I)/(l*l*l);
						K[1*n+2]=(2.0*E*I) / l;
			
						K[2*n+1]=-(6.0*E*I)/(l*l);
					
					}


					if(MPI_P_ID==0){
						disp(5,n, K, "K");
					}


					//SCALAPACK setup 

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
					desc_K[5]=5; //Local size
					desc_K[6]=0; //Reserved

					int desc_u_p[7]; 
					desc_u_p[0]=502; //Type is banded matrix P-by-1 (RHS)
					desc_u_p[1]=ctx; //Context
					desc_u_p[2]=3*(N_e-1); //Problem size
					desc_u_p[3]=3*((N_e-2)/2+1); //Blocking
					desc_u_p[4]=0; //Process row/column
					desc_u_p[5]=3*((N_e-2)/2+1); //Local size
					desc_u_p[6]=0; //Reserved
			
					int laf=(3*((N_e-2)/2+1)+2*4)*4; // laf
					int lwork=4*4; //lwork
					double *work=new double[lwork]; // work
					double *af=new double[laf]; // af
					double *af_tmp=new double[laf]; //af_tmp

					// Get and display the effective matrix [K_eff] ([K] is being overwritten with [K_eff])

					get_K_eff(dt,M,n,K);

					if(MPI_P_ID==0){
						disp(5,n,K,"K_eff");
					}

					//Convert K from row major to column major (for blas and lapack libraries)
					rm2cm(5,n,K);

					//Variables needed for solution
					double *u = new double[n](); //u at time layer n
					double *u_p = new double[n](); //u at time layer n+1
					double *u_t = new double[n](); //first derivative of u at time layer n
					double *u_t_p = new double[n](); //first derivative of u at time layer n+1
					double *u_tt = new double[n](); //second derivative of u at time layer n
					double *u_tt_p = new double[n](); //second derivative of u at time layer n+1
					double *tmp = new double[n](); //temporary storage used during computations
					int info=0;

					//Initialise time to 0
					t=0;

					//Display that we're solving implicit problem on two processes on process with ID 0
					if(MPI_P_ID==0){
						std::cout<<"Solving [M]d2{u}/dt2+[K]{u}={F} with implicit scheme on two processes"<<std::endl;
						std::cout<<std::endl;
					}

					F77NAME(pdpbtrf)('U',3*(N_e-1),4,K,1,desc_K,af,laf,work,lwork,info);       

					cblas_dcopy(laf,af,1,af_tmp,1); 

					if (info) {
						cout << "PDPBTRF error "<< info <<" !"<< endl;
					}  

					//Loop for Nt iterations
					for(int n_t=0;n_t<=N_t;n_t++){

						//Increment time
						t=t+dt;

						//Get F at iteration n+1
						get_F(t,T,l,F,N_n);

						//If process ID is 0 set last node y force to qy/l+Fy
						if(MPI_P_ID==0){
							F[3*(N_n-1)+1]=t*1000.0/T*(l+1);
						}

						//Calculate right hand side of equation
						for(int i=0;i<n;i++){
							tmp[i]=F[i]+M[i]*(1.0/(0.25*dt*dt)*u[i]+1.0/(0.25*dt)*u_t[i]+(1.0/(2.0*0.25)-1.0)*u_tt[i]);
						}

						//Copy the right hand side of the equation into u_p
						cblas_dcopy(n,tmp,1,u_p,1);
						
						//Solve the system [K]*{u_p}'={u_p} 

						F77NAME(pdpbtrs)('U',3*(N_e-1),4,1,K,1,desc_K,u_p,1,desc_u_p,af_tmp,laf,work,lwork,info);            

						if (info) {
							cout << "PDPBTRS error "<< info <<" !"<< endl;
						}  
						
						//Get u_tt_p and u_t_p
						for(int i=0;i<n;i++){
							u_tt_p[i]=1.0/(0.25*dt*dt)*(u_p[i]-u[i])-1.0/(0.25*dt)*u_t[i]-(1.0/(2.0*0.25)-1.0)*u_tt[i];  
							u_t_p[i]=u_t[i]+dt*0.5*u_tt[i]+dt*0.5*u_tt_p[i];
						}

						//Copy u_p into u, u_t_p into u_t and u_tt_p into u_tt
						cblas_dcopy(n,u_p,1,u,1);
						cblas_dcopy(n,u_t_p,1,u_t,1);
						cblas_dcopy(n,u_tt_p,1,u_tt,1);
						
					}

					Cblacs_gridexit(ctx);

					//Collect u results from processes
					MPI_Gather(u,n,MPI_DOUBLE,MPI_u,n,MPI_DOUBLE,0,MPI_COMM_WORLD);

					//Save collected u results to file and display them on process ID 0
					if(MPI_P_ID==0){
						
						N_n=N_e-1;
						n=3*N_n;

						log(N_n,l,MPI_u,"results/task_dynamic_implicit_parallel.log");

						cout<<endl<<"Done!"<<endl;

						disp(n,1,MPI_u,"u");

					}

				}

			}

		}

	//Finalise MPI
	MPI_Finalize();

	return 0;
}
