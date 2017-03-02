#include "functions.h"


void inv(double *A, int &r){


    int *IPIV = new int[r+1];
    int LWORK = r*r;
    double *WORK = new double[LWORK];
    int INFO;

    dgetrf_(&r,&r,A,&r,IPIV,&INFO);

    dgetri_(&r,A,&r,IPIV,WORK,&LWORK,&INFO);

    delete[] IPIV;
    delete[] WORK;

}
