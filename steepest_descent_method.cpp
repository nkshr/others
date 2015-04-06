#define ALFA 0.1
#include <iostream>

using namespace std;

double calc_func(double x)
{
  return x * x - 2 * x - 3;
}

double calc_grad(double x)
{
  return 2 * x - 2;
}

double update(double x)
{
  return x - ALFA * calc_grad(x);
}

bool check(double x)
{
  if(calc_grad(x) <= 0.01)
    return true;
  return false;
}
int main(int argc, char **argv)
{
  double x = 10;
  for(;;)
    {
      cout << calc_func(x) << endl;
      if(check(x))
	break;
      x = update(x);
    }
  cout << x << endl;
  return 0;
}
