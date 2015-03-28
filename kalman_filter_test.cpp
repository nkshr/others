#include <iostream>
#include <boost/random.hpp>
#include <ctime>
#include <memory.h>
#include <fstream>

using namespace std; 

const int num = 10;
int counter = 0;
double A[2][2] = { { 1, 0 }, { 0, 1 } };
double B[2][2] = { { 1, 0 }, { 0, 1 } };
double x[2] = { 0, 0 };
double z[2] = { 0, 0 };
double p_x[2] = { 0, 0 };
double p_z[2] = { 0, 0 };
double X[num][2];
double Z[num][2];
double U[num][2];
double u[2] = { 2, 2 };
double H[2][2] = { { 1, 0 }, { 0, 1 } };
double x_hat[2] = { 0, 0 };
double P[2][2] = {};
double Q[2][2] = { { 1, 0 }, { 0, 1 } };
double R[2][2] = { { 2, 0 }, { 0, 2 } };
double P_[2][2] = {};
double x_hat_[2] = {};
double temp0[2][2] = {};
double temp1[2][2] = {};
double K[2][2] = {};
double X_hat[num][2];

//generator
boost::random::mt19937 rng(static_cast<unsigned long>(time(0)));

//distributor
boost::random::normal_distribution<> nd0(0.0, 1.0);
boost::random::variate_generator< boost::random::mt19937,
	boost::random::normal_distribution<> > w(rng, nd0);

boost::random::normal_distribution<> nd1(0.0, 2.0);
boost::random::variate_generator< boost::random::mt19937,
	boost::random::normal_distribution<> > v(rng, nd1);


void get_current_state();
void get_measurement();
void time_update();
void measurement_update();
void apply_kf();
bool get_inverse(double X[2][2], double X_i[2][2]);
void show_matrix(double X[2][2]);

int main(int argc, char **argv)
{
	FILE *gp;
	gp = _popen("gnuplot -persist", "w");
	for (int i = 0; i < num; i++)
	{
		get_current_state();
		memcpy(X[i], x, sizeof(x));
		get_measurement();
		memcpy(Z[i], z, sizeof(z));
		memcpy(U[i], u, sizeof(u));
		apply_kf();
		memcpy(X_hat[i], x_hat, sizeof(x_hat));
	}

	ofstream t_ofs("true_data.txt");
	ofstream m_ofs("measurement_data.txt");
	ofstream e_ofs("estimate_data.txt");
	for (int i = 0; i < num; ++i)
	{
		t_ofs << X[i][0] << "	" << X[i][1] << endl;
		m_ofs << Z[i][0] << "	" << Z[i][1] << endl;
		e_ofs << X_hat[i][0] << "	" << X_hat[i][1] << endl;
	}
	fprintf(gp, "set xrange[-5:30]\n");
	fprintf(gp, "set yrange[-5:30]\n");
	fprintf(gp, "plot \"true_data.txt\" w lp pt 7\n");
	fprintf(gp, "replot \"measurement_data.txt\" w lp pt 7\n");
	fprintf(gp, "replot \"estimate_data.txt\" w lp pt 7\n");
	return 0;
}

void get_current_state()
{
	//x = A * x_ + B * u + w
	x[0] = A[0][0] * p_x[0] + A[0][1] * p_x[1] + 
		B[0][0] * u[0] + B[0][1] * u[1] + w();
	x[1] = A[1][0] * p_x[0] + A[1][1] * p_x[1] +
		B[1][0] * u[0] + B[1][1] * u[1] + w();
	memcpy(p_x, x, sizeof(x));
}

void get_measurement()
{
	//z = H * x + v
	z[0] = H[0][0] * x[0] + H[0][1] * x[1] + v();
	z[1] = H[1][0] * x[0] + H[1][1] * x[1] + v();
}

void apply_kf()
{
	time_update();
	measurement_update();
	counter++;
}

void time_update()
{
	//x^_ = A * x^ + B * u
	x_hat_[0] = A[0][0] * x_hat[0] + A[0][1] + x_hat[1] +
		B[0][0] * u[0] + B[0][1] + u[1];
	x_hat_[1] = A[1][0] * x_hat[0] + A[1][1] + x_hat[1] +
		B[1][0] * u[0] + B[1][1] + u[1];
	
	//P_ = A * P * A_t + Q 
	P_[0][0] = A[0][0] * (P[0][0] * A[0][0] + P[0][1] * A[0][1]) +
		A[0][1] * (P[1][0] * A[0][0] + P[1][1] * A[0][1]) + Q[0][0];
	P_[1][0] = A[1][0] * (P[0][0] * A[0][0] + P[0][1] * A[0][1]) +
		A[1][1] * (P[1][0] * A[0][0] + P[1][1] * A[0][1]) + Q[1][0];
	P_[0][1] = A[0][0] * (P[0][0] * A[1][0] + P[0][1] * A[1][1]) +
		A[0][1] * (P[1][0] * A[1][0] + P[1][1] * A[1][1]) + Q[0][1];
	P_[1][1] = A[1][0] * (P[0][0] * A[1][0] + P[0][1] * A[1][1]) +
		A[1][1] * (P[1][0] * A[1][0] + P[1][1] * A[1][1]) + Q[1][1];
}

void measurement_update()
{
	//temp0 = H * P_ * H_t + R
	temp0[0][0] = H[0][0] * (P_[0][0] * H[0][0] + P_[0][1] * H[0][1]) +
		H[0][1] * (P_[1][0] * H[0][0] + P_[1][1] * H[0][1]) + R[0][0];
	temp0[1][0] = H[1][0] * (P_[0][0] * H[0][0] + P_[0][1] * H[0][1]) +
		H[1][1] * (P_[1][0] * H[0][0] + P_[1][1] * H[0][1]) + R[1][0];
	temp0[0][1] = H[0][0] * (P_[0][0] * H[1][0] + P_[0][1] * H[1][1]) +
		H[0][1] * (P_[1][0] * H[1][0] + P_[1][1] * H[1][1]) + R[0][1];
	temp0[1][1] = H[1][0] * (P_[0][0] * H[1][0] + P_[0][1] * H[1][1]) +
		H[1][1] * (P_[1][0] * H[1][0] + P_[1][1] * H[1][1]) + R[1][1];
	show_matrix(temp0);
	if(!get_inverse(temp0, temp1))
		memcpy(temp1, temp0,sizeof(temp0));
	
	//K = P * H * temp1
	K[0][0] = P_[0][0] * (H[0][0] * temp1[0][0] + H[1][0] * temp1[1][0]) +
		P_[0][1] * (H[0][1] * temp1[0][0] + H[1][1] * temp1[1][0]);
	K[1][0] = P_[1][0] * (H[0][0] * temp1[0][0] + H[1][0] * temp1[1][0]) +
		P_[1][1] * (H[0][1] * temp1[0][0] + H[1][1] * temp1[1][0]);
	K[0][1] = P_[0][0] * (H[0][0] * temp1[0][1] + H[1][0] * temp1[1][1]) +
		P_[0][1] * (H[0][1] * temp1[0][1] + H[1][1] * temp1[1][1]);
	K[1][1] = P_[1][0] * (H[0][0] * temp1[0][1] + H[1][0] * temp1[1][1]) +
		P_[1][1] * (H[0][1] * temp1[0][1] + H[1][1] * temp1[1][1]);
	
	//x_hat = x_hat_ + K * (z - H * x_hat_)
	x_hat[0] = x_hat_[0] + K[0][0] * (Z[counter][0] - H[0][0] * x_hat_[0] - H[0][1] * x_hat_[1]) +
		K[0][1] * (Z[counter][1] - H[1][0] * x_hat_[0] - H[1][1] * x_hat_[1]);
	x_hat[1] = x_hat_[1] + K[1][0] * (Z[counter][0] - H[0][0] * x_hat_[0] - H[0][1] * x_hat_[1]) +
		K[1][1] * (Z[counter][1] - H[1][0] * x_hat_[0] - H[1][1] * x_hat_[1]);
	
	// P = P_ - K * H * P_
	P[0][0] = P_[0][0] -
		(K[0][0] * (H[0][0] * P_[0][0] + H[0][1] * P_[1][0]) +
		K[0][1] * (H[1][0] * P_[0][0] + H[1][1] * P_[1][0]));
	P[1][0] = P_[1][0] -
		(K[1][0] * (H[0][0] * P_[0][0] + H[0][1] * P_[1][0]) +
		K[1][1] * (H[1][0] * P_[0][0] + H[1][1] * P_[1][0]));
	P[0][1] = P_[0][1] -
		(K[0][0] * (H[0][0] * P_[0][1] + H[0][1] * P_[1][1]) +
		K[0][1] * (H[1][0] * P_[0][1] + H[1][1] * P_[1][1]));
	P[1][1] = P_[1][1] -
		(K[1][0] * (H[0][0] * P_[0][1] + H[0][1] * P_[1][1]) +
		K[1][1] * (H[1][0] * P_[0][1] + H[1][1] * P_[1][1]));
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

void show_matrix(double X[2][2])
{
	cout << X[0][0] << ", " << X[0][1] << endl;
	cout << X[1][0] << ", " << X[1][1] << endl;

}