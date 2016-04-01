// vp_1.cpp: определяет точку входа для консольного приложения.
//

#include "stdafx.h"
#include "Service.h"
#include "Jacobi.h"
#include <iostream>
#include <fstream>
#include<conio.h>

using namespace std;
const double a = 0.9;
const double b = 1.1;

bool isBorder(int i, int j, int N) {   // if knot of net is on border of my area
	if (i == 0 || j == 0 || j == N || (j <= N / 2 && i + j == N) || (j > N / 2 && i == j)) return true; else return false;
}

double f(double x, double y) {
	return (4 - 3.2 * x * x - 4.8 * y * y) / ((1 + x * x + y * y)*(1 + x * x + y * y)*(1 + x * x + y * y));
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
void solutionToMatrix(int ** Ind, double ** Result, double ** PreciseSolution ,double * x, int N, double step) { 
	/*fliis matrices Result and Precise solution with approximate solution and precise solution*/
	for (int i = 1; i < N + 1; i++) {
		for (int j = 1; j < N + 1; j++) {
			if (Ind[i][j] > 0) {
				PreciseSolution[i][j] = fi(step * i, step * j);
				Result[i][j] = x[Ind[i][j]];
			}
			else if (isBorder(i, j, N)) {
				Result[i][j] = fi(step * i, step * j);
				PreciseSolution[i][j] = fi(step * i, step * j);
			}
		}
	}
}

void fillMatrix(double ** A, int ** Ind, int n, double h, double * F, double step) { // Fills A matrix and F vector
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			if (Ind[i][j] > 0) {
				A[Ind[i][j]][Ind[i][j]] = 2 * (a + b) / h;
				F[Ind[i][j]] += f(step * i, step * j);
				if (Ind[i - 1][j] > 0) A[Ind[i][j]][Ind[i - 1][j]] = -a / h; else F[Ind[i][j]] += (a / h) * fi(step*(i - 1), step*j);
				if (Ind[i + 1][j] > 0) A[Ind[i][j]][Ind[i + 1][j]] = -a / h; else F[Ind[i][j]] += (a / h) * fi(step*(i + 1), step*j);
				if (Ind[i][j - 1] > 0) A[Ind[i][j]][Ind[i][j - 1]] = -b / h; else F[Ind[i][j]] += (b / h) * fi(step*i, step*(j - 1));
				if (Ind[i][j + 1] > 0) A[Ind[i][j]][Ind[i][j + 1]] = -b / h; else F[Ind[i][j]] += (b / h) * fi(step*i, step*(j + 1));

			}
		}
	}
}
template<typename T> T ** allocateMemory(int N) {
	T ** A;
	A = (T**)calloc(sizeof(T*), N);
	for (int i = 0; i <= N; i++)
		 A[i] = (T*)calloc(sizeof(T), N);
	return A;
}



int main()
{
	int N, k; // N is a grid size and k is quantity of variables in z vector
	cin >> N;
	double step = 1.0 / N;
	double h = step * step; 
	double **U, **Err;
	double *F, *x, *xNext;
	U = (double**)calloc(sizeof(double*), N + 1);
	for (int i = 0; i <= N; i++)
		U[i] = (double*)calloc(sizeof(double), N + 1);
	
	int **Ind;
	double **A, **Result, **PreciseSolution;

	Ind = allocateMemory<int>(N + 1);
	Result = allocateMemory<double>(N+1);
	PreciseSolution = allocateMemory<double>(N + 1);
	Err = allocateMemory<double>(N + 1);

	k = fillIndexes(Ind, N);

	cout << k << endl;
	F = (double *)calloc(sizeof(double), k + 1);
	x = (double *)calloc(sizeof(double), k + 1);
	xNext = (double *)calloc(sizeof(double), k + 1);
	x[1] = 1;
	A = allocateMemory<double>(k + 1); 
	

	fillMatrix(A, Ind,N, h, F, step);
	Solve(A, F, x, xNext, k);
	solutionToMatrix(Ind, Result, PreciseSolution, x, N, step);
	

	double max = normOfDiff(Result, PreciseSolution, Err, N);
	cout << max;

	print(Err, N, true, false, step);
	_getch();
	return 0;
}

