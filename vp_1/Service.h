#pragma once
#include <time.h>
#include <string>
using namespace std;

const double epsilon = 0.1;
struct Eigenvalues {
	double min, max;
};

void print(double * F, int n);


template <typename T> void print(T ** A, int n, bool toFile, bool withZero, double step) {
	if (toFile) {
		ofstream fout;
		fout.open("C:\\VP\\output.txt");
		for (int i = 0; i < n + 1; i++) {
			for (int j = 0; j < n + 1; j++) {
				fout << step * i << " " << step * j << " " << A[i][j] << endl;
			}
			fout << endl;
		}
	}
	else {
		for (int i = 0; i < n + 1; i++) {
			for (int j = 0; j < n + 1; j++) {
				cout << A[i][j] << " ";
			}
			cout << endl;
		}
		cout << endl;
	}

}
template<typename T> T ** allocateMemory(int N) {
	T ** A;
	A = (T**)calloc(sizeof(T*), N);
	for (int i = 0; i <= N; i++)
		A[i] = (T*)calloc(sizeof(T), N);
	return A;
}
void copy(double * from, double * to, int n);

double normOfDiff(double * x1, double * x2, int n);

double normOfDiff(double ** Result, double ** PreciseSolution, double ** Err, int N);

void mul(double ** aMatrix, double ** bMatrix, double ** product, int n);

void mul(double ** A, double * x, double * product, int n);

double norm(double * v, int n);