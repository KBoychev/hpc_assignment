///////////////////////////////////////////////////////////////////////////////////  
// rm2cm() Converts matrix from row major to column major
// ------------------------------------------------------------------------------
// @param r <int> - Rows
// @param c <int> - Columns
// @param m <double*> - Matrix/vector
// ------------------------------------------------------------------------------

#include "functions.h"


void rm2cm(int r,int c,double* m){


	double *m_tmp = new double[r*c]();
	
	for(int j=0;j<c;j++){
		for(int i=0;i<r;i++){
				m_tmp[j*r+i]=m[i*c+j];
		}
	}

	for(int k=0;k<r*c;k++){
		m[k]=m_tmp[k];
	}
	

	delete[] m_tmp;

}