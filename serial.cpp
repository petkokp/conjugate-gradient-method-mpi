#include <iostream>
#include <cmath>
#include <algorithm>
#include <stdio.h>
#include <stdlib.h>
#include <fstream>
#include <cmath>
#include <sys/time.h>
#include "shared.h"

using namespace std;
const double eps = 1e-6;

void matrixMulVec(double **a, double *b, double *r, int M, int N)
{
	for (int i = 0; i < M; i++)
	{
		r[i] = 0;
		for (int k = 0; k < N; k++)
		{
			r[i] += a[i][k] * b[k];
		}
	}
}

void solve(double **a, double *b, double *x, double *r, double *p, int M, int N)
{
	for (int i = 0; i < N; i++)
	{
		x[i] = 0;
	}
	double *tmp, *nr;
	tmp = new double[M];
	nr = new double[M];
	matrixMulVec(a, x, tmp, M, N);
	for (int i = 0; i < M; i++)
		r[i] = b[i] - tmp[i];
	for (int i = 0; i < M; i++)
		p[i] = r[i];
	int k = 0;
	while (true)
	{
		k++;
		matrixMulVec(a, p, tmp, M, N);
		double alpha = vecMulVec(r, p, M) / vecMulVec(p, tmp, M);
		for (int i = 0; i < N; i++)
			x[i] += alpha * p[i];
		for (int i = 0; i < M; i++)
			nr[i] = r[i] - alpha * tmp[i];
		if (check(nr, M))
			break;
		double beta = vecMulVec(nr, nr, M) / vecMulVec(r, r, M);
		for (int i = 0; i < M; i++)
			p[i] = nr[i] + beta * p[i];
		for (int i = 0; i < M; i++)
			r[i] = nr[i];
	}
	delete[] tmp;
	delete[] nr;
}

int main(int argc, char *argv[])
{
	int M, N;

	M = atoi(argv[1]);
	N = M;

	double elapsedTime, elapsedTime2;
	timeval start, end, end2;

	double **a, *b, *x, *r, *p;
	a = new double *[M];
	b = new double[M];
	x = new double[N];
	r = new double[M];
	p = new double[M];

	for (int i = 0; i < M; i++)
	{
		a[i] = new double[N];
	}

	ifstream matrixfile("matrix");
	if (!(matrixfile.is_open()))
	{
		cout << "Error: matrix file not found" << endl;
		return 0;
	}
	for (int i = 0; i < M; i++)
	{
		for (int j = 0; j < N; j++)
		{
			matrixfile >> a[i][j];
		}
	}
	for (int i = 0; i < M; i++)
	{
		matrixfile >> b[i];
	}
	matrixfile.close();

	gettimeofday(&start, NULL);
	solve(a, b, x, r, p, M, N);
	gettimeofday(&end, NULL);

	elapsedTime = (end.tv_sec - start.tv_sec) * 1000.0;
	elapsedTime += (end.tv_usec - start.tv_usec) / 1000.0;
	cout << "Time: " << elapsedTime << " ms." << endl
			 << endl;

	ofstream Af;
	Af.open("serial_results");
	for (int i = 0; i < N; i++)
	{
		Af << x[i] << endl;
	}
	Af.close();

	for (int i = 0; i < M; i++)
	{
		delete[] a[i];
	}
	delete[] a;
	delete[] b;
	delete[] x;
	delete[] r;
	delete[] p;
	return 0;
}