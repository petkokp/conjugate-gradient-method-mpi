#include <iostream>
#include <cmath>
#include <algorithm>
#include <stdio.h>
#include <stdlib.h>
#include <fstream>
#include <cmath>
#include <sys/time.h>
#include <mpi.h>
#include "shared.h"

using namespace std;
const double eps = 1e-6;
int myid, numprocs;
MPI_Status status;
int col;
double *t1, *t2, *res;
int k = 1;

int k1 = 100000, k2 = 200000, k3 = 300000;
int tt = 0;

void matrixMulVec(double **a, double *b, double *r, int M, int N)
{
	tt++;

	if (tt == 1)
	{
		for (int m = 1; m < numprocs; m++)
		{
			int size = 0;
			for (int j = (m - 1) * col; j < m * col; j++)
			{
				for (int i = 0; i < M; i++)
				{
					t1[size++] = a[i][j];
				}
			}

			MPI_Send(t1, size, MPI_DOUBLE, m, k1 + k * numprocs + m, MPI_COMM_WORLD);
		}
	}
	for (int m = 1; m < numprocs; m++)
	{
		int size2 = 0;
		for (int j = (m - 1) * col; j < m * col; j++)
		{
			t2[size2++] = b[j];
		}

		MPI_Send(t2, size2, MPI_DOUBLE, m, k2 + k * numprocs + m, MPI_COMM_WORLD);
	}

	for (int i = 0; i < M; i++)
	{
		r[i] = 0;
	}

	for (int m = 1; m < numprocs; m++)
	{
		MPI_Recv(res, M, MPI_DOUBLE, m, k3 + k * numprocs + m, MPI_COMM_WORLD, &status);
		for (int i = 0; i < M; i++)
		{
			r[i] += res[i];
		}
	}

	for (int j = (numprocs - 1) * col; j < N; j++)
	{
		for (int i = 0; i < M; i++)
		{
			r[i] += a[i][j] * b[j];
		}
	}
}

void Cal(int M, int N)
{
	tt++;
	if (tt == 1)
	{
		MPI_Recv(t1, col * M, MPI_DOUBLE, 0, k1 + k * numprocs + myid, MPI_COMM_WORLD, &status);
	}
	MPI_Recv(t2, col, MPI_DOUBLE, 0, k2 + k * numprocs + myid, MPI_COMM_WORLD, &status);
	for (int i = 0; i < M; i++)
	{
		res[i] = 0;
		for (int j = 0; j < col; j++)
		{
			res[i] += t1[j * M + i] * t2[j];
		}
	}

	MPI_Send(res, M, MPI_DOUBLE, 0, k3 + k * numprocs + myid, MPI_COMM_WORLD);
}

void solve(double **a, double *b, double *x, double *r, double *p, int M, int N)
{
	for (int i = 0; i < N; i++)
	{
		x[i] = rand() / (RAND_MAX + 1.0);
	}
	double *tmp, *nr;
	tmp = new double[M];
	nr = new double[M];
	matrixMulVec(a, x, tmp, M, N);
	for (int i = 0; i < M; i++)
		r[i] = b[i] - tmp[i];
	for (int i = 0; i < M; i++)
		p[i] = r[i];
	int flag;
	while (true)
	{
		k++;
		flag = 0;
		for (int m = 1; m < numprocs; m++)
		{
			MPI_Send(&flag, 1, MPI_INT, m, k * numprocs + m, MPI_COMM_WORLD);
		}
		matrixMulVec(a, p, tmp, M, N);
		double alpha = vecMulVec(r, r, M) / vecMulVec(p, tmp, M);
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
	flag = 1;
	k++;
	for (int m = 1; m < numprocs; m++)
	{
		MPI_Send(&flag, 1, MPI_INT, m, k * numprocs + m, MPI_COMM_WORLD);
	}
}

void solve2(double **a, double *b, double *x, double *r, double *p, int M, int N)
{
	int flag;
	Cal(M, N);
	while (true)
	{
		++k;
		MPI_Recv(&flag, 1, MPI_INT, 0, k * numprocs + myid, MPI_COMM_WORLD, &status);
		if (flag == 1)
			break;
		Cal(M, N);
	}
}

int main(int argc, char *argv[])
{

	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &myid);
	MPI_Comm_size(MPI_COMM_WORLD, &numprocs);

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

	col = N / numprocs;
	t1 = new double[M * col];
	t2 = new double[col];
	res = new double[M];

	if (myid == 0)
	{
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
		cout << elapsedTime << endl
				 << endl;

		ofstream Af;
		Af.open("parallel_results");
		for (int i = 0; i < N; i++)
		{
			Af << x[i] << endl;
		}
		Af.close();

		for (int i = 0; i < M; i++)
		{
			delete[] a[i];
		}
	}
	else
	{
		solve2(a, b, x, r, p, M, N);
	}
	delete[] a;
	delete[] b;
	delete[] x;
	delete[] r;
	delete[] p;
	delete[] t1;
	delete[] t2;
	delete[] res;
	MPI_Finalize();
	return 0;
}