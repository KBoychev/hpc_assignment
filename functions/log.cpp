#include "functions.h"


void log(int &n, double *u, std::string file) {

	std::ofstream log_file;
	log_file.open(file);

	if(log_file.good()){
		for(int i=0;i<n;i++){
			log_file<<u[i]<<"\n";
		}
	}

	log_file.close();

}