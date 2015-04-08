#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <memory.h>

using namespace std;

double x[7][2] = {{0.038, 0.050}, {0.194, 0.127},
		     {0.425, 0.094}, {0.626, 0.2122},
		     {1.253, 0.2729}, {2.500, 0.2665},
		     {3.740, 0.3317}};
int iter = 0;



double calc_gamma(double beta[2], double dst[7])
{
	for(int i = 0; i < 7; ++i)
		dst[i] = x[i][1] - ((beta[0]*x[i][0]) / (beta[1]+x[i][0]));
}

void calc_jacobian(double beta[2], double dst[7][2])
{
	for(int i = 0; i < 7; i++)
	{
		dst[i][0] = (-x[i][0]) / (beta[1]+x[i][0]);
		dst[i][1] = (beta[0]*x[i][0]) / pow(beta[1]+x[i][0], 2);
	}
}

void get_t_jacobian(double src[7][2], double dst[2][7])
{
	for(int i = 0; i < 7; i++)
	{
		dst[0][i] = src[i][0];
		dst[1][i] = src[i][1];
	}
}

bool get_inverse(double X[2][2], double X_i[2][2])
{
	if((X[0][0] * X[1][1] - X[0][1] * X[1][0]) < 0.1)
		return false;
	X_i[0][0] = (1 / (X[0][0] * X[1][1] - X[0][1] * X[1][0])) *
		X[1][1];
	X_i[0][1] = (1 / (X[0][0] * X[1][1] - X[0][1] * X[1][0])) *
		(-X[0][1]);
	X_i[1][0] = (1 / (X[0][0] * X[1][1] - X[0][1] * X[1][0])) *
		(-X[1][0]);
	X_i[1][1] = (1 / (X[0][0] * X[1][1] - X[0][1] * X[1][0])) *
		X[0][0];
	return true;
}

double update(double beta[2], double gamma[7], double jaco[7][2], double t_jaco[2][7], double dst[2])
{
	calc_gamma(beta, gamma);
	calc_jacobian(beta, jaco);
	get_t_jacobian(jaco, t_jaco);
	double temp[2][2]={}, temp_inv[2][2];
	for(int i = 0; i < 7; ++i)
	{
		temp[0][0] += t_jaco[0][i]*jaco[i][0];
		temp[1][0] += t_jaco[1][i]*jaco[i][0];
		temp[0][1] += t_jaco[0][i]*jaco[i][1];
		temp[1][1] += t_jaco[1][i]*jaco[i][1];
	}
	if(!get_inverse(temp, temp_inv))
		memcpy(temp_inv, temp, sizeof(temp));
	temp[0][0] = 0;
	temp[1][0] = 0;
	for(int i = 0; i < 7; ++i)
	{
		temp[0][0] += t_jaco[0][i]*gamma[i];
		temp[1][0] += t_jaco[1][i]*gamma[i];
	}
	dst[0] = beta[0] + temp_inv[0][0]*temp[0][0] + temp_inv[0][1]*temp[1][0];
	dst[1] = beta[1] + temp_inv[1][0]*temp[0][0] + temp_inv[1][1]*temp[1][0];
}

int main(int argc, char **argv)
{
	double jaco[7][2], t_jaco[2][7], gamma[7];
	double	 beta[2] = {0.9, 0.2}, dst[2] = {};
	FILE	*gp;
	gp		  = popen("gnuplot -persist", "w");
	fprintf(gp, "set title \"gauss newton method\"\n");
	fprintf(gp, "set xrange[0:4]\n");
	fprintf(gp, "set yrange[0:0.35]\n");
	int	 max_i = sizeof(x)/sizeof(x[0]);
	for(int i = 0; i < max_i; ++i)
	{
		fprintf(gp, "set label %d point pt 7 at %f,%f\n",
			i+1, x[i][0], x[i][1]);
	}
	for(int i = 0; i < 5; ++i)
	{
		update(beta, gamma, jaco, t_jaco, dst);
		memcpy(beta, dst, sizeof(dst));
		cout << beta[0] << " " << beta[1] << endl;
	}
	fprintf(gp, "f(x) = (%f*x) / (%f+x)\n", beta[0], beta[1]);
	fprintf(gp, "plot f(x)\n");
	return 0;
}
