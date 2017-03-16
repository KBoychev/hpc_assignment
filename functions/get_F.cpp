
#include "mpi.h"

#include "functions.h"


void get_F(double &t, double &T, double &l, double *F, int &N_n) {

	int MPI_N_P;
	int MPI_P_ID;

	MPI_Comm_size(MPI_COMM_WORLD,&MPI_N_P);
	MPI_Comm_rank(MPI_COMM_WORLD,&MPI_P_ID);


	for(int n_n=0;n_n<N_n;n_n++){

		F[3*n_n+0]=0;
		F[3*n_n+1]=t*1000.0/T*l;
		F[3*n_n+2]=0;

	}

	if(MPI_N_P==2 && MPI_P_ID==0){
		F[3*(N_n-1)+1]=t*1000.0/T*(l+1);
	}

	if(MPI_N_P==2 && MPI_P_ID==1){
		F[3*0+1]=t*1000.0/T*(l+1);
	}

	if(MPI_N_P==1){
		F[3*((N_n-1)/2)+1]=t*1000.0/T*(l+1);
	}

}
