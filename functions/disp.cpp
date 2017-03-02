#include "functions.h"


void disp(int r,int c, double *m , string l) {

	cout << endl;

	cout << setprecision(6);

	if (c == 1) {
		cout << r << "x" << c << " {" << l << "}=" << endl;
		for (int i = 0; i < r; i++) {
			cout << setw(12);
			cout << m[i];
			cout << endl;
		}
	}
	else {
		cout << r << "x" << c << " [" << l << "]=" << endl;

		for (int i = 0; i < r; i++) {
			for (int j = 0; j < c; j++) {
				cout << setw(12);
				cout << m[i*r + j];
			}
			cout << endl;
		}
	}
	cout << endl;
}