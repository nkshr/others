#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <iostream>

using namespace std;

#define W 10

double data[600][2];
double w[600];
double temp;

int get_random(int min, int max);
void apply_ls(double &a, double &b);
void e_estimate(double &a, double &b);

int main(int argc, char **argv)
{
	FILE *gp;
	srand((int)time(NULL));
	gp = _popen("gnuplot -persist", "w");
	fprintf(gp, "set title \"m_estimate\"\n");
	fprintf(gp, "set xrange[-10:330]\n");
	fprintf(gp, "set yrange[-10:330]\n");

	/* Make noise. */
	for (int i = 0; i < 300; i++)
	{
		data[i][0] = get_random(1, 300);
		data[i][1] = get_random(1, 300);
		while ((data[i][1] <= data[i][0] + 10 && 
			data[i][1] >= data[i][0] - 10))
		{
			data[i][1] = get_random(1, 300);
		}
	}

	/* Make true data. */
	for (int i = 0; i < 300; i++)
	{
		temp = (double)get_random(0, 10);
		if (i % 2 == 0)
			temp = (-1) * temp;
		data[i + 300][0] = i;
		data[i + 300][1] = i + temp;
	}

	double a, b;
	apply_ls(a, b);

	/* Draw a graph. */
	for (int i = 0; i < 300; i++)
	{
		fprintf(gp, "set label %d point pt 7 at %f,%f\n", i + 1, data[i][0], data[i][1]);
		fprintf(gp,
			"set label %d point pt 7 lc rgb \"spring-green\" at %f,%f\n",
			i + 301, data[i + 300][0], data[i + 300][1]);
	}
	fprintf(gp, "plot x\n");
	fprintf(gp, "replot %f * x + %f\n", a, b);

	e_estimate(a, b);
	fprintf(gp, "replot %f * x + %f\n", a, b);
	return 0;
}

int get_random(int min, int max)
{
	return min + (int)(rand()*(max - min + 1.0) / (1.0 + RAND_MAX));
}

void apply_ls(double &a, double &b)
{
	a = b = 0;
	//A = sigma (x * y), B = sigma (x), C = sigma (y), D = sigma (x ^ 2) 
	double A, B, C, D;
	A = B = C = D = 0;
	for (int i = 0; i < 600; ++i)
	{
		A += data[i][0] * data[i][1];
		B += data[i][0];
		C += data[i][1];
		D += pow(data[i][0], 2.0);
	}
	a = (600 * A - B * C) / (600 * D - pow(B, 2.0));
	b = (D * C - A * B) / (600 * D - pow(B, 2.0));
}

void e_estimate(double &a, double &b)
{
	for (int i = 0; i < 600; i++)
	{
		temp = data[i][1] - (a * data[i][0] + b);
		if (abs(temp) <= W)
			w[i] = pow(1 - pow((temp / W), 2), 2);
		else
			w[i] = 0;
	}

	a = b = 0;
	double A, B, C, D;
	A = B = C = D = 0;
	for (int i = 0; i < 600; ++i)
	{
		A += w[i] * (data[i][0] * data[i][1]);
		B += w[i] * data[i][0];
		C += w[i] * data[i][1];
		D += w[i] * (pow(data[i][0], 2.0));
	}
	a = (600 * A - B * C) / (600 * D - pow(B, 2.0));
	b = (D * C - A * B) / (600 * D - pow(B, 2.0));
}
