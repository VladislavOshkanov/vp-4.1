#include "stdafx.h"
#include "Lusternik.h"

#include <iostream>

Eigenvalues Lusternik_init(double ** A, int N) {
	Eigenvalues E;
	E.max = Lusternik(A, N);
	double ** B;
	B = allocateMemory<double>(N + 1);
	for (int i = 1; i < N + 1; i++)
		for (int j = 1; j < N + 1; j++) {
			B[i][j] = -A[i][j];
			if (i = j) B[i][j] += E.max;
		}
	E.min = E.max - Lusternik(B, N);

	return E;
}

double Lusternik(double ** A, int N) {
	double ** K, **L, **M;
	K = allocateMemory<double>(N + 1);
	L = allocateMemory<double>(N + 1);
	M = allocateMemory<double>(N + 1);
	double * y, *result;
	y = (double *)calloc(sizeof(double), N + 1);
	result = (double *)calloc(sizeof(double), N + 1);
	srand(time(NULL));
	for (int i = 1; i < N + 1; i++)
			y[i] = rand() % 10 - 5;
	copy(A, K, N);
	copy(A, M, N);
	double V, V_prev;
	V = 1;
	V_prev = 10;
	while (fabs(V - V_prev) > epsilon *  N) {
	//for (int q = 0; q < 10; q++){
		mul(K, M, L, N + 1);
		copy(L, K, N);

		mul(K, y, result, N + 1);
		double n1 = norm(result, N);
		mul(K, M, L, N + 1);
		mul(L, y, result, N + 1);
		double n2 = norm(result, N);
		V_prev = V;
		V = n2 / n1;
	}
	return V;
}