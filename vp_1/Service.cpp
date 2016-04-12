#include "stdafx.h"
#include<iostream>

#include "Service.h"

using namespace std;

double normOfDiff(double ** Result, double ** PreciseSolution, double ** Err, int N){ // norm of difference between two matrices
	double max = 0;
	for (int i = 1; i < N + 1; i++)
		for (int j = 1; j < N + 1; j++)
		{
			if (fabs(Result[i][j] - PreciseSolution[i][j]) > max) max = fabs(Result[i][j] - PreciseSolution[i][j]);
			Err[i][j] = fabs(Result[i][j] - PreciseSolution[i][j]);
		}
	return max;
}

void print(double * F, int n) {
	for (int i = 0; i < n; i++) {
		cout << F[i] << " ";
	}
	cout << endl;
}

void copy(double * from, double * to, int n) { // copies one vector to another
	for (int i = 1; i < n + 1; i++) {
		to[i] = from[i];
	}
}
void copy(double ** from, double ** to, int n) {
	for (int i = 1; i < n + 1; i++) {
		for (int j = 1; j < n + 1; j++) {
			to[i][j] = from[i][j];
		}
	}
}
double norm(double * v, int n) {
	double x = 0;
	for (int i = 1; i < n + 1; i++) {
		x += v[i] * v[i];
	}
	x = sqrt(x);
	return x;
}

void mul(double ** aMatrix, double ** bMatrix,double ** product, int n) {

	for (int row = 1; row < n; row++) {
		for (int col = 1; col < n; col++) {
			for (int inner = 1; inner < n; inner++) {
				product[row][col] += aMatrix[row][inner] * bMatrix[inner][col];
			}
		}
	}
}
void mul(double ** A, double * x, double * product, int n) {
	for (int row = 1; row < n; row++) {
		for (int inner = 1; inner < n; inner++) {
			product[row] += A[row][inner] * x[inner];
		}
	}
}
double normOfDiff(double * x1, double * x2, int n) { //norm of x1-x2 vectors
	double norm = 0;
	for (int i = 1; i < n + 1; i++) {
		norm += abs(x1[i] - x2[i]);
	}
	return norm;
}