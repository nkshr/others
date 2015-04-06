#include <iostream>

using namespace std;

double calc_deriv(double x)
{
  return 2 * x;
}

double calc_func(double x)
{
  return x * x - 2;
}

double update(double x)
{
  return x - (calc_func(x) / calc_deriv(x));
}

int main(int argc, char **argv)
{
  double x = 1000;
  for(;;)
    {
      cout << calc_func(x) << endl;
      if(calc_func(x) <= 0.001)
	break;
      x = update(x);
    }
  cout << x << endl;
  return 0;
}
