#include "stdafx.h"
#include "Jacobi.h"
#include "Service.h"

void Solve(double ** A, double * F, double * x, double * xNext, int N) { //solving matrix using Jacobi algorithm
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