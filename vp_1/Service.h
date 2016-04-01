#pragma once
#include <string>
using namespace std;

const double epsilon = 0.01;

void print(double * F, int n);


template <typename T> void print(T ** A, int n, bool toFile, bool withZero, double step) {
	if (toFile) {
		ofstream fout;
		fout.open("output.txt");
		for (int i = 0; i < n + 1; i++) {
			for (int j = 0; j < n; j++) {
				fout << step * i << " " << step * j << " " << A[i][j] << endl;
			}
			fout << endl;
		}
	}

}
void copy(double * from, double * to, int n);

double normOfDiff(double * x1, double * x2, int n);

double normOfDiff(double ** Result, double ** PreciseSolution, double ** Err, int N);