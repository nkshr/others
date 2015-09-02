#include "steepest_descent.h"
#include <string.h>
#include <stdlib.h>
#include <iostream>
using namespace std;

steepest_descent::steepest_descent(double  x[][2], int numx, double *y, double * init_params)
{
	this->x = (double **)malloc(sizeof(double)* numx);

	for (int i = 0; i < numx; ++i)
	{
		this->x[i] = (double*)malloc(sizeof(double)* 2);
	}

	alpha = init_params;
	
	for (int i = 0; i < numx; ++i)
	{
		memcpy(this->x[i], x[i], sizeof(double)* 2);
	}
	this->numx = numx;
	this->y = y;

	w = new double[2];
	w[0] = 0;
	w[1] = 0;

	gain =1;
}

void steepest_descent::update_params()
{
	double *temp = new double[numx];

	double gradient;

	double sum_gradient = 0;

	for (int i = 0; i < numx; ++i)
	{
		gradient = 1;

		for (int j = 0; j < numx; ++j)
		{
			gradient -= alpha[j] * y[j] * y[i] * (x[j][0] * x[i][0] + x[j][1] * x[i][1]);
		}

		sum_gradient += gradient * gradient;

		temp[i] = alpha[i] - gain * gradient;
	}

	memcpy(alpha, temp, sizeof(double)* numx);
	
	for (int i = 0; i < numx; ++i)
	{
		w[0] += alpha[i] * y[i] * x[i][0];
		w[1] += alpha[i] * y[i] * x[i][1];
	}

	cout << "gradient: " << gradient << endl;
	cout << w[0] << ", " << w[1] << endl;
	delete temp;
}

void steepest_descent::calc_Lagrangian()
{
	double sum_alpha = 0;

	for (int i = 0; i < numx; ++i)
	{
		sum_alpha += alpha[i];
	}

	double temp = 0;

	for (int i = 0; i < numx; ++i)
	{
		for (int j = 0; j < numx; ++j)
		{
			temp += alpha[i] * alpha[j] * y[i] * y[j] * (x[i][0] * x[j][0] + x[i][1] * x[j][1]);
		}
	}

	Lagrangian = sum_alpha + 0.5 * temp;
	cout << "Lagrangian: " << Lagrangian << endl;
}
