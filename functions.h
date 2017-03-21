

#ifndef _funcitons_h
#define _funcitons_h

	#include <string>
	
	////  
	// get_M() Populates the elements of the mass matrix [M]. 
	// ------------------------------------------------------------------------------
	// @param rho <double> - Density (kg/m^3)
	// @param A <double> - Cross-section area (m^2)
	// @param l <double> - Element length (m)
	// @param M <double*> - Mass matrix (Banded symmetric storage)
	// @param N_n <int> - Number of nodes
	// ------------------------------------------------------------------------------

	void get_M(double &rho, double &A, double &l, double *M, int &N_n);

	////  
	// get_K() Populates the elements of the stiffness matrix [K]. 
	// ------------------------------------------------------------------------------
	// @param E <double> - Young's modulus (Pa)
	// @param A <double> - Cross-section area (m^2)
	// @param I <double> - Second moment of area (m^4)
	// @param l <double> - Element length (m)
	// @param n <int> - Degrees of freedom
	// @param K <double*> - Stiffness matrix (Banded symmetric storage)
	// @param N_n <int> - Number of nodes
	// ------------------------------------------------------------------------------

	void get_K(double &E, double &A,  double &I, double &l, int &n, double *K, int &N_n);

	////  
	// get_F() Populates the elements of the force vector {F}. 
	// ------------------------------------------------------------------------------
	// @param t <double> - Time (s)
	// @param T <double> - Period (s)
	// @param l <double> - Element length (m)
	// @param F <double*> - Force vector 
	// @param N_n <int> - Number of nodes
	// ------------------------------------------------------------------------------

	void get_F(double &t, double &T, double &l, double *F, int &N_n);

	////  
	// get_K_eff() Populates the elements of the effective stiffness matrix [K_eff]. 
	// ------------------------------------------------------------------------------
	// @param dt <double> - Stepsize (s)
	// @param M <double*> - Mass matrix 
	// @param n <int> - Degrees of freedom
	// @param K <double*> - Stiffness matrix (Banded symmetric storage)
	// ------------------------------------------------------------------------------

	void get_K_eff(double &dt, double *M, int &n, double *K);


	////  
	// disp() Dispays matrix or vector. 
	// ------------------------------------------------------------------------------
	// @param r <int> - Rows
	// @param c <int> - Columns
	// @param m <double*> - Matrix/vector
	// @param l <string> - Matrix/vector label
	// ------------------------------------------------------------------------------

	void disp(int r,int c, double *m,std::string l);

	////  
	// log() Logs vertical displacments to file
	// ------------------------------------------------------------------------------
	// @param N_n <int> - Number of nodes
	// @param l <double> - Element length (m)
	// @param u <double*> - Displacements 
	// @param file <string> - File name
	// ------------------------------------------------------------------------------

	void log(int &N_n, double &l, double *u, std::string file);

	////  
	// rmjr() Converts matrix from row major to column major
	// ------------------------------------------------------------------------------
	// @param r <int> - Rows
	// @param c <int> - Columns
	// @param m <double*> - Matrix/vector
	// ------------------------------------------------------------------------------

	void rmjr(int r,int c,double *m);

#endif
