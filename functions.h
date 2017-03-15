

#ifndef _funcitons_h
#define _funcitons_h

	#include <iostream>
	#include <iomanip>
	#include <string>
	#include <fstream>

	void get_M(double &rho, double &A, double &l, double *M, int &N_n);

	void get_K(double &E, double &A,  double &I, double &l, int &n, double *K, int &N_n);

	void get_F(double &t, double &T, double &l, double *F, int &N_n);

	void get_K_eff(double &dt, double *M, int &n, double *K);

	void disp(int,int,double*,std::string);

	void log(int &n, double *u, std::string file);

#endif
