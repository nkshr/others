#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

using namespace std;

double x[7][2] = {{0.038, 0.050}, {0.194, 0.127},
		     {0.425, 0.094}, {0.626, 0.2122},
		     {1.253, 0.2729}, {2.500, 0.2665},
		     {3.740, 0.3317}};
int iter = 0;

double update(double x, double y)
{
  
}

double calc_gamma(double b1, double b2, double dst[7])
{
	for(int i = 0; i < 7; ++i)
		dst[i] = x[i][1] - ((b1*x[i][0]) / (b2+x[i][0]));
}

void calc_jacobian(double b1, double b2, double dst[7][2])
{
	for(int i = 0; i < 7; i++)
	{
		dst[i][0] = (-x[i][0]) / (b2+x[i][0]);
		dst[i][1] = pow((b1*x[i][0]) / (b2+x[i][0]), 2);
	}
}

void get_t_jacobian(double src[7][2], double dst[2][7])
{
	for(int i = 0; i < 7; i++)
	{
		dst[0][i] = src[i][0];
		dst[0][i] = src[i][1];
	}
}

int main(int argc, char **argv)
{
	double jaco[7][2], t_jaco[2][7], gamma[7];
	double	 beta1	  = 1, beta2 = 1;
	FILE	*gp;
	gp		  = popen("gnuplot -persist", "w");
	fprintf(gp, "set title \"gauss newton method\"\n");
	fprintf(gp, "set xrange[0:4]\n");
	fprintf(gp, "set yrange[0:0.35]\n");
	int	 max_i = sizeof(x)/sizeof(x[0]);
	for(int i = 0; iter < max_i; ++i)
	{
		fprintf(gp, "set label %d point pt 7 at %f,%f\n",
			i+1, x[i][0], x[i][1]);
	}
	calc_gamma(beta1, beta2, gamma);
	calc_jacobian(beta1, beta2, jaco);
	get_t_jacobian(jaco, t_jaco);
	fprintf(gp, "f(x) = (%f*x) / (%f+x)\n", beta1, beta2);
	fprintf(gp, "plot f(x)\n");
	return 0;
}
