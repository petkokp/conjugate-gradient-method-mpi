#include "shared.h"

static double eps = 1e-6;
int t = 0;

bool check(double *r, int N)
{
    double res = 0;
    for (int i = 0; i < N; i++)
    {
        res += std::fabs(r[i]);
    }
    t++;

    if (res < eps)
        return true;
    else
        return false;
}

double vecMulVec(double *a, double *b, int N)
{
    double ans = 0;
    for (int i = 0; i < N; i++)
    {
        ans += a[i] * b[i];
    }
    return ans;
}
