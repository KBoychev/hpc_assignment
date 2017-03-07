

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
  
 
 
  	int N_e;
	double l;
	double T;
	int N_t;
  	double dt;
  	int r_p=3;
  	double *M_p=new double[r_p*r_p]();
  	double *K_p=new double[r_p*r_p]();
  	double *K_pp=new double[r_p*r_p]();
  	double *K_pm=new double[r_p*r_p]();
  	double *F_p=new double[r_p]();
  	double *u_p;
  	double *u_pp;
  	double *u_pm;
  	double *u_buff=new double[r_p]();

  	double t;
  	

 	if(MPI_P_ID==0){
		
		double L = 10.0; 
		double A = 0.1*0.12; 
		double I = (0.1*0.12*0.12*0.12)/12; 
		double E = 210000000000; 
		double rho = 7850; 


		N_e = 24; 	
		l =L/N_e;
		T = 1; 
		N_t = 10000;
		dt=T/(N_t-1);


		u_p=new double[(N_e/2+2)*3]();
		u_pp=new double[(N_e/2+2)*3]();
		u_pm=new double[(N_e/2+2)*3]();

		MPI_Send(&N_e,1,MPI_INT,1,0,MPI_COMM_WORLD);	
		MPI_Send(&l,1,MPI_DOUBLE,1,0,MPI_COMM_WORLD);	
		MPI_Send(&T,1,MPI_DOUBLE,1,0,MPI_COMM_WORLD);	
		MPI_Send(&N_t,1,MPI_INT,1,0,MPI_COMM_WORLD);	
		MPI_Send(&dt,1,MPI_DOUBLE,1,0,MPI_COMM_WORLD);	

		int r_e=6;

		double *M_e = new double[r_e*r_e]();

		double *K_e = new double[r_e*r_e]();

		set_M_e(r_e,r_e, M_e, rho, A, l);

		disp(r_e,r_e,M_e,"M_e");

		set_K_e(r_e, r_e,K_e, E, A, I, l);

		disp(r_e,r_e,K_e,"K_e");

		
		M_p[0*r_p+0]=2*M_e[0*r_e+0];
		M_p[1*r_p+1]=2*M_e[1*r_e+1];
		M_p[2*r_p+2]=2*M_e[2*r_e+2];

		//Cool stuff happens here. We don't need the large stiffness and mass
		//matrices! The displacements at each node are dependent on the stiffness matrix of the
		//node and the stiffness matrix of the adjacent nodes, therefore we can use the element 
		//matrix only and extract the node stiffnesses matrices from there :D no memory needed for large matrices!

		MPI_Send(M_p,r_p*r_p,MPI_DOUBLE,1,0,MPI_COMM_WORLD);

		K_p[0*r_p+0]=K_e[0*r_e+0]+K_e[3*r_e+3];
		K_p[1*r_p+1]=K_e[1*r_e+1]+K_e[4*r_e+4];
		K_p[2*r_p+2]=K_e[2*r_e+2]+K_e[5*r_e+5];

		MPI_Send(K_p,r_p*r_p,MPI_DOUBLE,1,0,MPI_COMM_WORLD);

		K_pp[0*r_p+0]=K_e[0*r_e+3];
		K_pp[1*r_p+1]=K_e[1*r_e+4];
		K_pp[1*r_p+2]=K_e[1*r_e+5];
		K_pp[2*r_p+1]=K_e[2*r_e+4];
		K_pp[2*r_p+2]=K_e[2*r_e+5];

		MPI_Send(K_pp,r_p*r_p,MPI_DOUBLE,1,0,MPI_COMM_WORLD);

		K_pm[0*r_p+0]=K_e[3*r_e+0];
		K_pm[1*r_p+1]=K_e[4*r_e+1];
		K_pm[1*r_p+2]=K_e[4*r_e+2];
		K_pm[2*r_p+1]=K_e[5*r_e+1];
		K_pm[2*r_p+2]=K_e[5*r_e+2];
		
		MPI_Send(K_pm,r_p*r_p,MPI_DOUBLE,1,0,MPI_COMM_WORLD);

		delete[] K_e;

		delete[] M_e;

	}


	//// Set the process matrices M_p, K_p
	//--------------------------------------------------------

	if(MPI_P_ID!=0){

  		MPI_Recv(&N_e,1,MPI_DOUBLE,0,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);

  		u_p=new double[(N_e/2+2)*3]();
		u_pp=new double[(N_e/2+2)*3]();
		u_pm=new double[(N_e/2+2)*3]();


  		MPI_Recv(&l,1,MPI_DOUBLE,0,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
		MPI_Recv(&T,1,MPI_DOUBLE,0,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
		MPI_Recv(&N_t,1,MPI_INT,0,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
		MPI_Recv(&dt,1,MPI_DOUBLE,0,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
		MPI_Recv(M_p,r_p*r_p,MPI_DOUBLE,0,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
		MPI_Recv(K_p,r_p*r_p,MPI_DOUBLE,0,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
		MPI_Recv(K_pp,r_p*r_p,MPI_DOUBLE,0,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
		MPI_Recv(K_pm,r_p*r_p,MPI_DOUBLE,0,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
		
	}

	//// Solve the dynamic problem [M]d2{u}/dt2+[K]{u}={F}. Solution
	//	 is performed with explicit integration method.
	//--------------------------------------------------------
	

	if(MPI_P_ID==0){
		std::cout<<"Solving [M]d2{u}/dt2+[K]{u}={F} on "<<MPI_NP<<" processes."<<std::endl;
		std::cout<<std::endl;
	}

	

	for(int n_t=1;n_t<=N_t;n_t++){

		if(MPI_P_ID==0){

			std::cout<<"Iteration "<<n_t<<" t="<<t<<std::endl;

		}

		F_p[0]=0;
		F_p[2]=0;


		if(MPI_P_ID==0){

			for(int k=1;k<=(N_e/2);k++){

				if(k==(N_e/2)){
					F_p[1]=t*1000.0/T*l+t*1000.0/T;
				}else{
					F_p[1]=t*1000.0/T*l;
				}

				for(int i=0;i<3;i++){

					u_pp[k*3+i]=0;

					for(int j=0;j<3;j++){

						u_pp[k*3+i]=u_pp[k*3+i]-K_pm[i*r_p+j]*u_p[(k-1)*3+j]-K_p[i*r_p+j]*u_p[k*3+j]-K_pp[i*r_p+j]*u_p[(k+1)*3+j];
					}

					u_pp[k*3+i]=u_pp[k*3+i]+F_p[i]+2.0/(dt*dt)*M_p[i*r_p+i]*u_p[k*3+i]-1.0/(dt*dt)*M_p[i*r_p+i]*u_pm[k*3+i];
					u_pp[k*3+i]=dt*dt*1.0/M_p[i*r_p+i]*u_pp[k*3+i];

				}
			}

			//send to proceess 1

			u_buff[0]=u_pp[(N_e/2-1)*3+0];
			u_buff[1]=u_pp[(N_e/2-1)*3+1];
			u_buff[2]=u_pp[(N_e/2-1)*3+2];


			MPI_Send(u_buff,3,MPI_DOUBLE,1,0,MPI_COMM_WORLD);
			
		}

		if(MPI_P_ID==1){

			//receive from proceess 0

			MPI_Recv(u_buff,3,MPI_DOUBLE,0,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);

			u_pp[0]=u_buff[0];
			u_pp[1]=u_buff[1];
			u_pp[2]=u_buff[2];


			for(int k=1;k<=(N_e/2);k++){

				if(k==1){
					F_p[1]=t*1000.0/T*l+t*1000.0/T;
				}else{
					F_p[1]=t*1000.0/T*l;
				}

				for(int i=0;i<3;i++){

					u_pp[k*3+i]=0;

					for(int j=0;j<3;j++){

						u_pp[k*3+i]=u_pp[k*3+i]-K_pm[i*r_p+j]*u_p[(k-1)*3+j]-K_p[i*r_p+j]*u_p[k*3+j]-K_pp[i*r_p+j]*u_p[(k+1)*3+j];
					}

					u_pp[k*3+i]=u_pp[k*3+i]+F_p[i]+2.0/(dt*dt)*M_p[i*r_p+i]*u_p[k*3+i]-1.0/(dt*dt)*M_p[i*r_p+i]*u_pm[k*3+i];
					u_pp[k*3+i]=dt*dt*1.0/M_p[i*r_p+i]*u_pp[k*3+i];

				}
			}

			//send to process 0

			u_buff[0]=u_pp[6];
			u_buff[1]=u_pp[7];
			u_buff[2]=u_pp[8];

			MPI_Send(u_buff,3,MPI_DOUBLE,0,0,MPI_COMM_WORLD);

		}


		if(MPI_P_ID==0){

			//receive from process 1

			MPI_Recv(u_buff,3,MPI_DOUBLE,1,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);

			u_pp[(N_e/2+1)*3+0]=u_buff[0];
			u_pp[(N_e/2+1)*3+1]=u_buff[1];
			u_pp[(N_e/2+1)*3+2]=u_buff[2];

		
		}



		for(int i=0;i<((N_e/2+2)*3);i++){
			u_pm[i]=u_p[i];
			u_p[i]=u_pp[i];
		}

		t=t+dt;
	}

	//Gather/reduce results
	//-------------------------------

	if(MPI_P_ID==0){

		std::cout<<"Done!"<<std::endl;
		std::cout<<std::endl;

		disp((N_e/2+2)*3,1,u_pp,"u_pp");

	}

	MPI_Finalize();

	return 0;
}
