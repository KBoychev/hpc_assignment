

#ifndef _funcitons_h
#define _funcitons_h

	#include <iostream>
	#include <iomanip>
	#include <string>
	#include <fstream>

	
	
	// void set_M_e(int&,int&,double*,double&,double&,double&); //? do we need this

	// void set_K_e(int&,int&,double*,double&,double&,double&,double&); //? do we need this

	// void get_M(int&,int&,double*,int&,double*,int&);

	// void get_K(int&,int&,double*,int&,double*, int&);

	// void get_K_eff(int&, double*, double*, double*, double&); 

	// void get_F(int&, double*, int&, double&, double&);


	void get_M(double &rho, double &A, double &l, double *M, int &N_n);

	void get_K(double &E, double &A,  double &I, double &l, int &n_k, double *K, int &N_n);

	void get_F(double &t, double &T, double &l, double *F, int &N_n);

	void get_K_eff(double &dt, double *M, int n_k, double *K);


	// void inv(int&,double*);

	void disp(int,int,double*,std::string);

	void log(int &n, double *u, std::string file);

#endif
