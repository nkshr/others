#pragma comment(lib, "opencv_core2410d")

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <time.h>
#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>

using namespace std;
using namespace cv;

double data[900][2];
double inliers[300][2];
int num_data = sizeof(data) / sizeof(data[0]);
int num_inliers = sizeof(inliers) / sizeof(inliers[0]);
int it = 0;
double threshold = 8000000;

int get_random(int min, int max);
void solve_equations(double &a, double &b);
void re_estimate(double &a, double &b);
void push_inlier(double x, double y);

int main(int argc, char **argv)
{
	FILE *gp;
	srand((int)time(NULL));
	gp = _popen("gnuplot -persist", "w");
	//fprintf(gp, "set size square\n");
	double temp;
	fprintf(gp, "set title \"ransac test\"\n");
	fprintf(gp, "set xrange[-10:330]\n");
	fprintf(gp, "set yrange[-10:330]\n");

	/* Make noise. */
	for (int i = 0; i < 300; i++)
	{
		data[i][0] = get_random(1, 300);
		data[i][1] = get_random(1, 300);
		while ((data[i][1] <= data[i][0] + 10 && data[i][1] >= data[i][0] - 10)
			|| (data[i][1] <= (0.5 * data[i][0]) + 10 && data[i][1] >= (0.5 * data[i][0]) - 10))
		{
			data[i][1] = get_random(1, 300);
		}
	}

	/* Make true data1. */
	for (int i = 0; i < 300; i++)
	{
		temp = (double)get_random(0, 10);
		if (i % 2 == 0)
			temp = (-1) * temp;
		data[i + 300][0] = i;
		data[i + 300][1] = i + temp; 
	}

	/* Make true data2. */
	for (int i = 0; i < 300; i++)
	{
		temp = (double)get_random(0, 10);
		if (i % 2 == 0)
			temp = (-1) * temp;
		data[i + 600][0] = i;
		data[i + 600][1] = 0.5 * i + temp;
	}

	/* Draw a graph. */
	for (int i = 0; i < 300; i++)
	{
		fprintf(gp, "set label %d point pt 7 at %f,%f\n", i + 1, data[i][0], data[i][1]);
		fprintf(gp, 
			"set label %d point pt 7 lc rgb \"spring-green\" at %f,%f\n", 
			i + 301, data[i + 300][0], data[i + 300][1]);
		fprintf(gp,
			"set label %d point pt 7 lc rgb \"spring-green\" at %f,%f\n",
			i + 601, data[i + 600][0], data[i + 600][1]);
	}
	fprintf(gp, "plot x, 0.5 * x\n");

	/* Solve simiultaniou equations.: 
		a * x1 + b = y1
		a * x2 + b = y2 */
	double a1, b1, a2, b2;
	for (int line_counter = 0; line_counter < 2;)
	{
		if (line_counter == 0)
		{
			solve_equations(a1, b1);
			line_counter++;
			continue;
		}
		solve_equations(a2, b2);
		if (fabs(a2 - a1) > 0.4)
			line_counter++;
	}
	fprintf(gp, "replot %f * x + %f\n", a1, b1);
	fprintf(gp, "replot %f * x + %f\n", a2, b2);
	fprintf(gp, "set terminal png\n");
	fprintf(gp, "set out \"ransac.png\"\n");
	fprintf(gp, "replot");
	return 0;
}

int get_random(int min, int max)
{
	return min + (int)(rand()*(max - min + 1.0) / (1.0 + RAND_MAX));
}

void solve_equations(double &a, double &b)
{
	Mat left_side, right_side, solution;
	int x1, x2, counter;
	double p_range;	//a permissible range
	p_range = 30;
	bool found = false;
	it = 0;
	for (int i = 0; i < 373; i++)
	{
		x1 = get_random(0, 899);
		x2 = get_random(0, 899);
		solution = Mat(2, 1, CV_64FC1);
		left_side = (Mat_<double>(2, 2) << data[x1][0], 1, data[x2][0], 1);
		right_side = (Mat_<double>(2, 1) << data[x1][1], data[x2][1]);
		solve(left_side, right_side, solution);
		a = solution.at<double>(0, 0);
		b = solution.at<double>(1, 0);
		counter = 0;
		for (int j = 0; j < num_data; j++)
		{
			if (data[j][1] <= (a * data[j][0] + b + 10) && data[j][1] >= (a * data[j][0] + b - 10))
			{
				push_inlier(data[j][0], data[j][1]);
				counter++;
			}
		}
		if (counter > num_inliers)
		{
			re_estimate(a, b);
			break;
		}
	}
}

void re_estimate(double &a, double &b)
{
	a = b = 0;
	//A = sigma (x * y), B = sigma (x), C = sigma (y), D = sigma (x ^ 2) 
	double A, B, C, D;
	A = B = C = D = 0;
	for (int i = 0; i < num_inliers; ++i)
	{
		A += inliers[i][0] * inliers[i][1];
		B += inliers[i][0];
		C += inliers[i][1];
		D += pow(inliers[i][0], 2.0);
	}
	a = (num_inliers * A - B * C) / (num_inliers * D - pow(B, 2.0));
	b = (D * C - A * B) / (num_inliers * D - pow(B, 2.0));
}

void push_inlier(double x, double y)
{
	inliers[it][0] = x;
	inliers[it][1] = y;
	++it;
}