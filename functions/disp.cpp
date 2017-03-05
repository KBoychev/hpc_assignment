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

		for (int i = 0; i < r; i++) {
			for (int j = 0; j < c; j++) {
				std::cout << std::setw(12);
				std::cout << m[i*c + j];
			}
			std::cout << std::endl;
		}
	}
	std::cout << std::endl;
}