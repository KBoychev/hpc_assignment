////  
// log() Logs vertical displacments to file
// ------------------------------------------------------------------------------
// @param N_n <int> - Number of nodes
// @param l <double> - Element length (m)
// @param u <double*> - Displacements 
// @param file <string> - File name
// ------------------------------------------------------------------------------

#include <fstream>

#include "functions.h"


void log(int &N_n, double &l, double *u, std::string file) {

	std::ofstream log_file;
	log_file.open(file);

	if(log_file.good()){
		for(int n_n=0;n_n<N_n;n_n++){
			log_file<<(n_n+1)*l<<","<<u[3*n_n+1]<<"\n";
		}
	}

	log_file.close();

}