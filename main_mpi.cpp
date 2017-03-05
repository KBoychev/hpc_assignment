

#include <iostream>
#include <iomanip>




#include "mpi.h"
#include "functions.h"

int main(int argc, char* argv[]) {

	int MPI_NP;
	int MPI_P_ID;
	int MPI_status;


	MPI_status = MPI_Init(&argc, &argv);

	MPI_status = MPI_Comm_size(MPI_COMM_WORLD,&MPI_NP);

	if(MPI_NP!=2){
		return 1;
	}

  	MPI_status = MPI_Comm_rank(MPI_COMM_WORLD,&MPI_P_ID);
  
  	//std::cout<<"Process "<<MPI_P_ID+1<<" of "<<MPI_NP<<std::endl;

  	int N_e_p;
  	int N_e;
	double l;
	double T;
	int N_t;
  	double dt;
  	double *M;
  	double *K;
  	double t;
  	int r;
  	int r_e;
  	double qy;
  	double Fy;

  	
 	if(MPI_P_ID==0){
		
		double L = 10.0; 
		double A = 0.1*0.12; 
		double I = (0.1*0.12*0.12*0.12)/12; 
		double E = 210000000000; 
		double rho = 7850; 


		N_e = 4; 
		N_e_p=N_e/2;
		l =L/N_e;
		T = 1; 
		N_t = 10000;
		dt=T/(N_t-1);

		MPI_Send(&N_e,1,MPI_INT,1,0,MPI_COMM_WORLD);	
		MPI_Send(&l,1,MPI_DOUBLE,1,0,MPI_COMM_WORLD);	
		MPI_Send(&T,1,MPI_DOUBLE,1,0,MPI_COMM_WORLD);	
		MPI_Send(&N_t,1,MPI_INT,1,0,MPI_COMM_WORLD);	
		MPI_Send(&dt,1,MPI_DOUBLE,1,0,MPI_COMM_WORLD);	

		r_e = 6; 

		double *M_e = new double[r_e*r_e]();

		double *K_e = new double[r_e*r_e]();

		set_M_e(r_e,r_e, M_e, rho, A, l);

		set_K_e(r_e, r_e,K_e, E, A, I, l);

		r = 3 * (N_e_p + 1); 

		K=new double[r*r]();
	
		M= new double[r*r]();
	
		get_M(r, r, M, r_e, M_e, N_e_p);

		//disp(r,r,M,"M");

		MPI_Send(M,r*r,MPI_DOUBLE,1,0,MPI_COMM_WORLD);

		get_K(r, r, K, r_e, K_e, N_e_p);

		MPI_Send(K,r*r,MPI_DOUBLE,1,0,MPI_COMM_WORLD);

		//disp(r,r,K,"K");

		delete[] K_e;

		delete[] M_e;

	}


	//// Set the process matrices M_p, K_p
	//--------------------------------------------------------

	if(MPI_P_ID!=0){


  		MPI_Recv(&N_e,1,MPI_DOUBLE,0,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
  		MPI_Recv(&l,1,MPI_DOUBLE,0,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
		MPI_Recv(&T,1,MPI_DOUBLE,0,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
		MPI_Recv(&N_t,1,MPI_INT,0,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
		MPI_Recv(&dt,1,MPI_DOUBLE,0,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);

		N_e_p=N_e/2;
		r_e=6;
		r=3 * (N_e_p + 1);
  		

  		K=new double[r*r]();
  		M=new double[r*r]();


		MPI_Recv(M,r*r,MPI_DOUBLE,0,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);

		MPI_Recv(K,r*r,MPI_DOUBLE,0,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);

		
	}

	//// Solve the dynamic problem [M]d2{u}/dt2+[K]{u}={F}. Solution
	//	 is performed with explicit integration method.
	//--------------------------------------------------------
	


	if(MPI_P_ID==0){
		//std::cout<<"Solving [M]d2{u}/dt2+[K]{u}={F}"<<std::endl;
		//std::cout<<std::endl;
	}

	

	
	

	double *F_e = new double[r_e];
	double *F = new double[r];
	//double *F_p = new double[r]();
	double *u= new double[r*3](); 
	double *tmp = new double[r]();
	double *buff = new double[3]();



	for(int n_t=1;n_t<=N_t;n_t++){

		if(MPI_P_ID==0){

			std::cout<<"Iteration "<<n_t<<" t="<<t<<std::endl;

		}

		//Get F (needs optimisation)
		//-------------------------------

		qy=t*1000.0/T;
		Fy=0;

		set_F_e(r_e,F_e,l,qy);

		get_F(r,F,r_e,F_e,N_e_p,Fy);

		Fy=t*1000.0/T;

		if(MPI_P_ID==0){
			F[0]=0;
			F[1]=0;
			F[2]=0;
			F[r-2]=F[r-2]+Fy/2.0;
		}

		if(MPI_P_ID==1){
			F[1]=F[1]+Fy/2.0;
			F[r-3]=0;
			F[r-2]=0;
			F[r-1]=0;
		}

		

		//Get u at next timestep
		//-------------------------------

		for(int i=0;i<r;i++){

			tmp[i]=0;

			for(int j=0;j<r;j++){

				tmp[i]=tmp[i]-(K[i*r+j]-2.0/(dt*dt)*M[i*r+j])*u[j*3+1]-1.0/(dt*dt)*M[i*r+j]*u[j*3+2];

			}

			tmp[i]=tmp[i]+F[i];

			u[i*3]=tmp[i]*dt*dt*1.0/M[i*r+i];
		}

		

		//Send overlapping u
		//-------------------------------

		//Send u[r-3] from process 0 to process 1
		//--------------------------------------------------
		if(MPI_P_ID==0){

			buff[0]=u[(r-3)*3];
			buff[1]=u[(r-2)*3];
			buff[2]=u[(r-1)*3];

			MPI_Send(buff,3,MPI_DOUBLE,1,0,MPI_COMM_WORLD);

			//std::cout<<"Process 0 sent data!"<<std::endl;

			//disp(3,1,buff,"buff");
		}

		//Receive u[r-3] u[r-2] u[r-1] from process 0 and send u[2] to process 0
		//--------------------------------------------------
		if(MPI_P_ID==1){

			MPI_Recv(buff,3,MPI_DOUBLE,0,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);

			//std::cout<<"Process 1 received data!"<<std::endl;

			u[0*3]=u[0*3]+buff[0];
			u[1*3]=u[0*3]+buff[1];
			u[2*3]=-u[0*3];//+buff[2];

			buff[0]=u[0*3];
			buff[1]=u[1*3];
			buff[2]=u[2*3];


			MPI_Send(buff,3,MPI_DOUBLE,0,0,MPI_COMM_WORLD);

			//std::cout<<"Process 1 sent data!"<<std::endl;
		}

		//Receive u[2] from process 1
		//--------------------------------------------------
		if(MPI_P_ID==0){

			MPI_Recv(buff,3,MPI_DOUBLE,1,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);

			u[(r-3)*3]=buff[0];
			u[(r-2)*3]=buff[1];
			u[(r-1)*3]=buff[2];

			//std::cout<<"Process 0 received data!"<<std::endl;

		}

		//Set boundary conditions 
	    //-------------------------------

		if(MPI_P_ID==0){
			u[0*3]=0;
			u[1*3]=0;
			u[2*3]=0;
			//u[(r-1)*3]=0;
		}

		if(MPI_P_ID==1){
			//u[2*3]=0;
			u[(r-1)*3]=0;
			u[(r-2)*3]=0;
			u[(r-3)*3]=0;
		}

		//Save u values 
	    //-------------------------------

		for(int i=0;i<r;i++){
			u[i*3+2]=u[i*3+1];
			u[i*3+1]=u[i*3];
		}

		

		t=t+dt;
	}


	//Gather/reduce results
	//-------------------------------

	if(MPI_P_ID==0){

		std::cout<<"Done!"<<std::endl;
		std::cout<<std::endl;
		disp(r,3,u,"u");

	}


	MPI_Finalize();

	return 0;
}
