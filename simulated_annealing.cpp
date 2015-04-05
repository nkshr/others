#define ALFA 0.5

#include <iostream>

using namespace std;

double annealing(double startState, int maxIter, double goalE);
double EVAL(double state);
double NEIGHBOUR(double state);
int get_random(int min, int max);
double PROBABILITY(double e1, double e2, double t);
double TEMPERATURE(double r);
double random();

int main(int argc, char **argv)
{
	FILE *gp;
	gp = _popen("gnuplot -persist", "w");
	fprintf(gp, "set title \"simulated annealing\"\n");
	fprintf(gp, "set xrange[-4:8]\n");
	fprintf(gp, "set yrange[-50:15]\n");
	fprintf(gp, "set grid\n");
	fprintf(gp, "plot x*x*x*x - 7*x*x*x + 5*x*x + 31*x - 30\n");
	double result;
	for (;;)
	{
		result = annealing(100, 1000000, -400);
		cout << "result: " << result << endl;
		break;
	}
	return 0;
}

double annealing(double startState, int maxIter, double goalE)
{
	double state, bestState, nextState, e, bestE, nextE;
	state = startState;
	e = EVAL(state);
	bestState = state;
	bestE = e;
	for (int iter = 0; iter < maxIter; iter++)
	{
		nextState = NEIGHBOUR(state);
		nextE = EVAL(nextState);
		if (nextE < bestE)
		{
			bestState = nextState;
			bestE = nextE;
			if (bestE <= goalE)
				return bestState;
		}
		if (random() <= PROBABILITY(e, nextE, TEMPERATURE((double)iter / maxIter)))
		{
			state = nextState;
			e = nextE;
		}
	}
	return bestState;
}

double EVAL(double state)
{
	double x = state;
	return pow(x, 4) - 7 * pow(x, 3) + 5 * pow(x, 2) + 31 * x - 30;
}

double NEIGHBOUR(double state)
{
	if (get_random(-10, 9) >= 0)
		return state + 0.1;
	else
		return state - 0.1;
}

int get_random(int min, int max)
{
	return min + (int)(rand()*(max - min + 1.0) / (1.0 + RAND_MAX));
}

double PROBABILITY(double e1, double e2, double t)
{
	//cout << t << endl;
	if (e1 >= e2)
		return 1;
	else
		return exp((e1 - e2) / t);
}

double random()
{
	return 0.1 * get_random(0, 1);
}

double TEMPERATURE(double r)
{
	return pow(ALFA, r);
}