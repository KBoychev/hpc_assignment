

#include "mpi.h"

#include "functions.h"


void get_K(double &E, double &A,  double &I, double &l, int &n, double *K, int &N_n) {


	int MPI_N_P;
	int MPI_P_ID;

	MPI_Comm_size(MPI_COMM_WORLD,&MPI_N_P);
	MPI_Comm_rank(MPI_COMM_WORLD,&MPI_P_ID);

	for(int n_n=0;n_n<N_n;n_n++){	

		if(n_n<N_n-1){

			K[0*n+(3*n_n+5)]=(6.0*E*I)/(l*l);
			K[1*n+(3*n_n+3)]=-(A*E)/l;
			K[1*n+(3*n_n+4)]=-(12.0*E*I)/(l*l*l);
			K[1*n+(3*n_n+5)]=(2.0*E*I) / l;
			K[2*n+(3*n_n+4)]=-(6.0*E*I)/(l*l);		
		}


		K[4*n+(3*n_n)]=2*(A*E)/l;
		K[4*n+(3*n_n+1)]=2*(12.0*E*I)/(l*l*l);
		K[4*n+(3*n_n+2)]=2*(4.0*E*I)/l;
	}

	if(MPI_N_P==2 && MPI_P_ID==0){
		K[4*n+(3*(N_n-1)+0)]=(A*E)/l;
		K[4*n+(3*(N_n-1)+1)]=(12.0*E*I)/(l*l*l);
		K[4*n+(3*(N_n-1)+2)]=(4.0*E*I)/l;
	}
	
	if(MPI_N_P==2 && MPI_P_ID==1){
		K[4*n+0]=(A*E)/l;
		K[4*n+1]=(12.0*E*I)/(l*l*l);
		K[4*n+2]=(4.0*E*I)/l;
	}

}