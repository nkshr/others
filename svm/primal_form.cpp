#include <iostream>
#include <fstream>
#include <string.h>
#include <stdlib.h>
using namespace std;

#include <time.h>

#include "steepest_descent.h"

int get_random(int main, int max);

double X1[30][2], X2[30][2];


int main(int argc, char ** argv)
{
	srand((int)time(NULL));
	FILE *gp;

#ifdef WIN32
	gp = _popen("gnuplot -persist", "w");
#else
	gp = popen("gnuplot -persist", "w");
#endif

	ofstream X1_ofs("./X1.txt");
	ofstream X2_ofs("./X2.txt");

	bool detected1 = false, detected2 = false;

	for (int i = 0; i < 30; ++i)
	{
		while (!detected1)
		{
			X1[i][0] = (double)get_random(0, 9) + (double)get_random(0, 9) * 0.1;
			X1[i][1] = (double)get_random(0, 9) + (double)get_random(0, 9) * 0.1;

			if (X1[i][1] > X1[i][0]+1.414)
				detected1 = true;
		}

		while (!detected2)
		{
			X2[i][0] = (double)get_random(0, 9) + (double)get_random(0, 9) * 0.1;
			X2[i][1] = (double)get_random(0, 9) + (double)get_random(0, 9) * 0.1;

			if (X2[i][1] < X2[i][0] - 1.414)
				detected2 = true;
		}
		
		X1_ofs << X1[i][0] << "\t" << X1[i][1] << endl;
		X2_ofs << X2[i][0] << "\t" << X2[i][1] << endl;

		detected1 = false;
		detected2 = false;
	}

	double X[60][2], Y[60], alpha[60];

	for (int i = 0; i < 30; ++i)
	{
		Y[i] = 1;
		Y[i + 30] = -1;
	}

	for (int i = 0; i < 60; ++i)
		alpha[i] = 1;

	memcpy(X, X1, sizeof(X1[0]) * 30);
	memcpy(X+30, X2, sizeof(X2[0]) * 30);
	steepest_descent sd(X, 60, Y, alpha);

	for (int i = 0; i < 20; ++i)
	{
		sd.update_params();
		sd.calc_Lagrangian();
	}

	fprintf(gp, "set xrange[0:10]\n");
	fprintf(gp, "set yrange[0:10]\n");
	fprintf(gp, "plot \"X1.txt\"  pt 7\n");
	fprintf(gp, "replot \"X2.txt\"  pt 7\n");

	return 0;
}

int get_random(int min, int max)
{
	return min + (int)(rand()*(max - min + 1.0) / (1.0 + RAND_MAX));
}


