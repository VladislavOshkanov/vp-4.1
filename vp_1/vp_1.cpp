// vp_1.cpp: определяет точку входа для консольного приложения.
//

#include "stdafx.h"
#include <iostream>
#include <fstream>
#include<conio.h>
using namespace std;
const double a = 0.9;
const double b = 1.1;
const double epsilon = 0.01;
double f(double x, double y) {
	return (4 - 3.2 * x * x - 4.8 * y * y) / ((1 + x * x + y * y)*(1 + x * x + y * y)*(1 + x * x + y * y));
}
template<typename T>void print(T ** A, int n, bool toFile, bool withZero){
	int start;
	if (withZero) start = 0; else start = 1;
	if (!toFile) {
		for (int i = start; i < n; i++) {
			for (int j = start; j < n; j++)
				cout << A[i][j] << " ";
			cout << endl;
		}
	}
}
void print(double * F, int n) {
	for (int i = 0; i < n; i++) {
		cout << F[i] << " ";
	}
	cout << endl;
}
double normOfDiff (double * x1, double * x2, int n) {
	double norm = 0;
	for (int i = 1; i < n + 1; i++) {
		norm += abs(x1[i] - x2[i]);
	}
	return norm;
}
void copy(double * from, double * to, int n) {
	for (int i = 1; i < n + 1; i++) {
		to[i] = from[i];
	}
}
void Solve(double ** A, double * F, double * x, double * xNext, int N) {
	double sum = 0;
	while (normOfDiff(x, xNext, N) > epsilon) {
		copy(xNext, x, N);
		for (int i = 1; i < N + 1; i++) {
			sum = 0;
			for (int j = 1; j < N + 1; j++) {
				if (i != j) sum += A[i][j] * x[j];
			}
			xNext[i] = (1 / A[i][i])*(F[i] - sum);
		}
	}
}
template<typename T> void setZero(T ** A, int n) {
	for (int i = 0; i < n; i++)
		for (int j = 0; j < n; j++)
			A[i][j] = 0;
}
int fillIndexes(int ** Ind, int N) {
	int k = 0;
	for (int i = 0; i <= N; i++)
	{
		for (int j = 0; j <= N; j++)
		{
			if ((i > 0 && i < N && j > 0 && j < N) && ((j <= N / 2 && i + j < N) || (j > N / 2 && i < j))) {
				k++;
				Ind[i][j] = k;
			}
		}
	}
	return k;

}
double fi(double x, double y) {
	return 1 / (1 + x*x + y*y);

}
bool isBorder(int i, int j, int N) {
	if (i == 0 || j == 0 || j == N || (j <= N / 2 && i + j == N) || (j > N / 2 && i == j)) return true; else return false;
}
void fillMatrix(double ** A, int ** Ind, int k, int n, double h, double * F, double step) {
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			if (Ind[i][j] > 0) {
				A[Ind[i][j]][Ind[i][j]] = 2 * (a + b) / h;
				F[Ind[i][j]] += f(step * i, step * j);
				//F[Ind[i][j]] = Ind[i][j];
				if (Ind[i - 1][j] > 0) A[Ind[i][j]][Ind[i - 1][j]] = -a / h; else F[Ind[i][j]] += (a / h) * fi(step*(i - 1), step*j); // TODO add f(x,y) value
				if (Ind[i + 1][j] > 0) A[Ind[i][j]][Ind[i + 1][j]] = -a / h; else F[Ind[i][j]] += (a / h) * fi(step*(i + 1), step*j);
				if (Ind[i][j - 1] > 0) A[Ind[i][j]][Ind[i][j - 1]] = -b / h; else F[Ind[i][j]] += (b / h) * fi(step*i, step*(j - 1));
				if (Ind[i][j + 1] > 0) A[Ind[i][j]][Ind[i][j + 1]] = -b / h; else F[Ind[i][j]] += (b / h) * fi(step*i, step*(j + 1));

			}
		}
	}
}



int main()
{
	int N, k; // N is a grid size and k is quantity of variables in z vector
	cin >> N;
	double step = 1.0 / N;
	double h = step * step; 
	double **U;
	double *F, *x, *xNext;
	U = (double**)calloc(sizeof(double*), N + 1);
	for (int i = 0; i <= N; i++)
		U[i] = (double*)calloc(sizeof(double), N + 1);
	
	int **Ind;
	double **A;
	Ind = (int**)calloc(sizeof(int*), N + 1);
	for (int i = 0; i <= N; i++)
		Ind[i] = (int*)calloc(sizeof(int), N + 1);
	setZero(Ind, N + 1);
	k = fillIndexes(Ind, N);
	print(Ind, N + 1, false, true);
	cout << k << endl;
	F = (double *)calloc(sizeof(double), k + 1);
	x = (double *)calloc(sizeof(double), k + 1);
	xNext = (double *)calloc(sizeof(double), k + 1);
	x[1] = 1;
	A = (double**)calloc(sizeof(double*), k + 1);
	for (int i = 0; i <= k + 1; i++)
		A[i] = (double*)calloc(sizeof(double), k + 1);
	fillMatrix(A, Ind, k, N, h, F, step);
	Solve(A, F, x, xNext, k);
	print(x, k);
	//print(F, k);
	cout << fi(step * 2, step * 2) << endl;
	//print(A, k + 1, false, false);
	//cout << (normOfDiff(x, xNext, k)) << endl;
	_getch();
	return 0;
}

