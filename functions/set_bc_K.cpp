#include "functions.h"


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