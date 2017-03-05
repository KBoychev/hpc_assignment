

#ifndef _funcitons_h
#define _funcitons_h

	#include <iostream>
	#include <iomanip>
	#include <string>

	#include <Accelerate/Accelerate.h>

	void set_M_e(int&,int&,double*,double&,double&,double&);

	void set_K_e(int&,int&,double*,double&,double&,double&,double&);

	void set_F_e(int&,double*,double&, double&);

	void get_M(int&,int&,double*,int&,double*,int&);

	void get_M_p(int, int&,int&,double*,double*);

	void get_K(int&,int&,double*,int&,double*, int&);

	void get_K_p(int,int&,int&,double*,double*);

	void get_F(int&,double*,int&,double*,int&,double&);

	void get_F_p(int,int&,double*,int&,double*);

	void inv(double*,int&);

	void disp(int,int,double*,std::string);

#endif
