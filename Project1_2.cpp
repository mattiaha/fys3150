#include <cstdlib>
#include <iostream>
#include <cmath>
#include <iomanip>
#include "time.h" 
#include "armadillo"
using namespace std;
using namespace arma;
double f(double x) {
	return 100.0 * exp(-10.0 * x);
}
double exact(double x) { return 1.0 - (1 - exp(-10)) * x - exp(-10 * x); }
int main(int argc, char* argv[])

{
	clock_t start, finish;
	int n;
	cout << "dimension of nxn matrix" << endl;
	cin >> n;
	
	double h = 1.0 / (n + 1.0);
	double hh = h * h;
	vec a(n), b(n), c(n),x(n),v(n),bi(n);
	mat AO = randu<mat>(n, n);
	for (int i = 0; i < n; i++) {
		x(i) = h * i;
		bi(i) = hh * f(x(i));
		for (int j = 0; j < n; j++) if (i != j) {
			if (i - j != 1) {
				if (j - i != 1) {

					AO(i, j) = 0;
				}
			}
		}
	}

	start = clock();
	a(0) = 0;
	b(0) = AO(0,0);
	c(n - 1) = 0;
	
	for (int o = 1; o < n; o++) {
		a(o) = AO(o,o-1);
		b(o) = AO(o, o);
		c(o-1) = AO(o-1,o);
	}
	

	vec cdot(n), ddot(n);
	cdot(0) = c(0) / b(0);
	ddot(0) = bi(0) / b(0);
	for (int s=1; s< n; s++){
		cdot(s) = c(s)/(b(s)-cdot(s-1)*a(s));
		ddot(s) = (bi(s)-ddot(s-1)*a(s))/(b(s)-cdot(s-1)*a(s));
	}
	v(n - 1) = ddot(n - 1);
	for (int sa = n-2; sa > 0; sa--) {
		v(sa) = ddot(n - 1) - cdot(sa) * v(sa + 1);
	}
	finish = clock();
	vec cf(n);
	for (int p = 0; p < n; p++) {
		cf(p) = exact(x(p));
	}
	v(0) = 0;
	double timeused = (finish - start) / ((double)CLOCKS_PER_SEC);
	cout << setprecision(10) << setw(20) << "Time used for computation=" << timeused << endl;
	// removing the comment marks in the following section allows the values to be written into a text file, and then be plotted in python.
	//ofstream myfile("1000.txt");
	//cout << v << endl;
	//for (int count = 0; count < n; count++) {
		//myfile << v(count) << " ";
	//}
	//for (int coun = 0; coun < n; coun ++){
		//myfile << cf(coun) << " ";
	//}
	//myfile.close();
	return 0;


}