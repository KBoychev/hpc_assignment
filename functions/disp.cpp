
#include <iostream>
#include <iomanip>
#include <cmath>

#include "functions.h"


void disp(int r,int c, double *m , std::string l) {

	std::cout << std::endl;

	std::cout << std::setprecision(6);

	if (c == 1) {
		std::cout << r << "x" << c << " {" << l << "}=" << std::endl;
		for (int i = 0; i < r; i++) {
			std::cout << std::setw(12);
			std::cout << m[i];
			std::cout << std::endl;
		}
	}
	else {
		
		std::cout << r << "x" << c << " [" << l << "]=" << std::endl;
		std::cout << std::endl;

		int k_min;
		int k_max;

		for(int k=0;k<=floor(c/10);k++){

			k_min=10*k;
			k_max=(k+1)*10;

			if(k==0){
				k_min=0;
				k_max=10;
			}

			if(k==floor(c/10)){
				k_min=10*k;
				k_max=c;
			}

			std::cout<<"Columns "<<k_min+1<<" through "<<k_max<<std::endl;
			std::cout << std::endl;

			for (int i = 0; i < r; i++) {

				for (int j = k_min; j < k_max; j++) {
					std::cout << std::setw(12);
					std::cout << m[i*c+j];
				}

				std::cout << std::endl;
			}

			std::cout << std::endl;
		}
			

		
	}
	std::cout << std::endl;
}