#include "stdafx.h"
#include "Lusternik.h"
#include"Service.h"

double Lusternik_init(double ** A, int N) {

	return 0;
}

double Lusternik(double ** A, int N) {
	double ** K, **L, **M;
	K = allocateMemory<double>(N + 1);
	L = allocateMemory<double>(N + 1);
	M = allocateMemory<double>(N + 1);
	double * y, *result;
	y = (double *)calloc(sizeof(double), 4);
	result = (double *)calloc(sizeof(double), 4);
	srand(time(NULL));
	/*for (int i = 1; i < 4; i++)
		for (int j = 1; j < 4; j++) {
			K[i][j] = rand() % 10 - 5;
			M[i][j] = K[i][j];
			y[j] = rand() % 10 - 5;
		}*/
	for (int i = 1; i < 40; i++) {
		mul(K, M, L, 5);
		for (int j = 1; j < 5; j++)
			for (int l = 1; l < 5; l++)
				K[j][l] = L[j][l];
	}
	mul(K, y, result, 4);
	double n1 = norm(result, 4);
	mul(K, M, L, 4);
	mul(L, y, result, 4);
	double n2 = norm(result, 4);

	return 0;
}