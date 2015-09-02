class steepest_descent
{
private:
	int numx;
	//double *(x[2]);
	double **x;
	double *w, *y, *alpha, *params;
	double Lagrangian, gain;
public:
	steepest_descent(double  x[][2], int num_xd, double *y, double * init_params);
	void update_params();
	void calc_Lagrangian();
};