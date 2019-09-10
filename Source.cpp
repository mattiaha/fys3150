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
	double hh = log(h) * log(h);
	vec bi(n), x(n), v(n);
	for (int i = 0; i < n; i++) {
		x(i) = h * i;
		bi(i) = hh * f(x(i));
	}


	start = clock();
	vec c(n), d(n);
	c(0) = -1.0 / 2.0;
	d(0) = bi(0) / 2.0;
	for (int j = 1; j < n; j++) {
		c(j) = -1.0 / (2 - c(j - 1));
		d(j) = (bi(j) + d(j - 1) / (2.0 + c(j - 1)));
	}
	v(n - 1) = d(n - 1);
	for (int k = n - 2; k > -1; k--) {
		v(k) = d(k) - c(k) * v(k + 1);
	}
	finish = clock();
	vec cf(n);
	for (int p = 0; p < n; p++) {
		cf(p) = exact(x(p));
	}
	vec r_error(n);
	r_error(0) = 0;
	for (int e = 1; e < n; e++) {
		r_error(e) = (log10(abs((v(e) - cf(e)) / cf(e))));
	}


	double timeused = (finish - start) / ((double)CLOCKS_PER_SEC);
	cout << setprecision(10) << setw(20) << "Time used for computation=" << timeused << endl;
	ofstream myfile("10000000.txt");

	for (int coun = 0; coun<n;coun ++) {
		myfile << r_error(coun) << " ";
	}
	myfile.close();



	return 0;

}