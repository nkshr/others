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

int main(int argc, char ** argv)
{
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
	double square_error;
	for (int j = 0; j < 10; ++j)
	{
		calc_jacobian(jacobian, x_vec, t_vec);
		//cout << jacobian << endl;
		for (int i = 0; i < y_vec.rows; ++i)
		{
			temp1.at<double>(i, 0) = y_vec.at<double>(i, 0)
				- calc_func(x_vec, t_vec.at<double>(i, 0));
		}
		temp0 = jacobian.t() * jacobian;
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
		
		cout << delta << endl;
		cout << x_vec << endl;
		cout << jacobian << endl;
		cout << square_error << endl << endl;
		if (square_error < 10)
			break;
	}
	fprintf(gp, "set xrange[0:10]\n");
	fprintf(gp, "set yrange[0:20]\n");
	fprintf(gp, "plot \"data.txt\" pt 7\n");
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