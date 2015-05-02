#include <iostream>
#include <opencv2/opencv.hpp>
#include <opencv2/opencv_lib.hpp>
#include <string>
#include <fstream>

using namespace std;
using namespace cv;

double calc_func(Mat &x_vec, double t);
void calc_jacobian(Mat &jacobian, Mat &x_vec, Mat &t_vec);
int get_random(int min, int max);
double get_metric(double square_error, double pre_square_error,
	Mat delta, double rambda, Mat jacobian, Mat difference);

int main(int argc, char ** argv)
{
	ofstream ofs("square_error.txt");
	Mat y_vec(5, 1, CV_64F), t_vec(5, 1, CV_64F), x_vec(2, 1, CV_64F), 
		delta(2, 1, CV_64F), jacobian(5, 2, CV_64F);
	ifstream ifs("data.txt");
	string str;
	FILE *gp;
	gp = _popen("gnuplot -persist", "w");
	if (ifs.fail() == true)
	{
		cerr << "Couldn't open the fail." << endl;
		exit(1);
	}
	for (int i = 0; i < 5; ++i)
	{
		getline(ifs, str);
		t_vec.at<double>(i, 0) = stod(str.substr(0, (int)str.rfind(" ")));
		y_vec.at<double>(i, 0) = stod(str.substr((int)str.rfind(" ") + 1, -1));
	}
	x_vec.at<double>(0, 0) = 1;
	x_vec.at<double>(1, 0) = 1;
	Mat temp0, temp1(5, 1, CV_64F);
	double square_error, pre_square_error = 0, rambda = 0.1, metric, 
		accept_threshold = 0.5, L_d = 9, L_u = 11;
	for (int j = 0; j < 20; ++j)
	{
		calc_jacobian(jacobian, x_vec, t_vec);
		//cout << jacobian << endl;
		for (int i = 0; i < y_vec.rows; ++i)
		{
			temp1.at<double>(i, 0) = y_vec.at<double>(i, 0)
				- calc_func(x_vec, t_vec.at<double>(i, 0));
		}
		temp0 = jacobian.t() * jacobian + rambda * jacobian.t() * jacobian;
		if (determinant(temp0) == 0)
		{
			delta.at<double>(0, 0) = (double)get_random(-10, 10) * 0.1;
			delta.at<double>(1, 0) = (double)get_random(-10, 10) * 0.1;

		}
		else
		{
			delta = temp0.inv() * jacobian.t() * temp1;
		}
		x_vec += delta;

		/*Calculate square error.*/
		for (int i = 0; i < t_vec.rows; i++)
		{
			square_error = 0;
			square_error +=
				pow((y_vec.at<double>(i, 0) -
				calc_func(x_vec, t_vec.at<double>(i, 0))), 2.0);
		}
		ofs << j << " " << square_error << endl;
		metric = get_metric(square_error, pre_square_error, delta, rambda, jacobian, temp1);
		if (j == 0)
		{
			pre_square_error = square_error;
			continue;
		}
		cout << metric << endl;

		if (metric > accept_threshold)
		{
			(rambda / L_d) > pow(10, -7.0) ? rambda = rambda / L_d : rambda = pow(10, -7);
			pre_square_error = square_error;
			cout << "accept" << " " << j << endl;
		}
		else
		{
			x_vec -= delta;
			(rambda * L_u) < pow(10, 7) ? rambda = rambda * L_u : rambda = pow(10, 7);
		}
	}
	fprintf(gp, "set xrange[0:10]\n");
	fprintf(gp, "set yrange[0:20]\n");
	fprintf(gp, "plot \"data.txt\" pt 7\n");

	FILE *gp2 = _popen("gnuplot -persist", "w");
	fprintf(gp2, "plot \"square_error.txt\" w lp pt 7\n");
	fprintf(gp, "replot %f * exp(%f * x)\n", 2.5, 0.25);
	fprintf(gp, "replot %f * exp(%f * x)\n", x_vec.at<double>(0, 0), x_vec.at<double>(1, 0));
	/*cout << t_vec << endl;
	cout << y_vec << endl;*/
	return 0;
}

double calc_func(Mat &x_vec, double t)
{
	return  x_vec.at<double>(0, 0) * exp(x_vec.at<double>(1, 0) * t);
}

double calc_func(double x1, double x2, double t)
{
	return  x1 * exp(x2 * t);
}

void calc_jacobian(Mat &jacobian, Mat &x_vec, Mat &t_vec)
{
	for (int i = 0; i < t_vec.rows; ++i)
	{
		/*jacobian.at<double>(i, 0) = 
			calc_func((double)(x_vec.at<double>(0, 0) + 0.01), x_vec.at<double>(1, 0), t_vec.at<double>(i, 0))
			- calc_func(x_vec.at<double>(0, 0), x_vec.at<double>(1, 0), t_vec.at<double>(i, 0));
		jacobian.at<double>(i, 0) /= 0.01;
		jacobian.at<double>(i, 1) =
			calc_func(x_vec.at<double>(0, 0), (double)(x_vec.at<double>(1, 0) + 0.01), t_vec.at<double>(i, 0))
			- calc_func(x_vec.at<double>(0, 0), x_vec.at<double>(1, 0), t_vec.at<double>(i, 0));
		jacobian.at<double>(i, 1) /= 0.01;*/
		jacobian.at<double>(i, 0) = exp(x_vec.at<double>(1, 0) * t_vec.at<double>(i, 0));
		jacobian.at<double>(i, 1) = x_vec.at<double>(0, 0) * t_vec.at<double>(i, 0) * 
			exp(x_vec.at<double>(1, 0) * t_vec.at<double>(i, 0));
	}
}

int get_random(int min, int max)
{
	return min + (int)(rand()*(max - min + 1.0) / (1.0 + RAND_MAX));
}

double get_metric(double square_error, double pre_square_error, 
	Mat delta, double rambda, Mat jacobian, Mat difference)
{
	Mat metric =
	(pre_square_error - square_error)
		/ (2 * delta.t() * (rambda * delta + jacobian.t() * difference));
	return metric.at<double>(0, 0);
}