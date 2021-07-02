#include <cmath>
#include <stdexcept>
#include "w.h"

using namespace std;

double w_2(double x)
{
  auto abs_x = abs(x);
  if (abs_x >= 0 && abs_x <= 1)
  {
    return 1 - abs_x;
  }
  return 0;
}

double w_3(double x)
{
  auto abs_x = abs(x);

  if (abs_x >= 0 && abs_x <= 1)
  {
    return 1.0 - ((5.0 / 2.0) * pow(abs_x, 2.0)) + ((3.0 / 2.0) * pow(abs_x, 3.0));
  }
  else if (abs_x >= 1 && abs_x <= 2)
  {
    return ((1.0 / 2.0) * pow(2.0 - abs_x, 2.0)) * (1.0 - abs_x);
  }

  return 0;
}

double w_4(double x)
{
  auto abs_x = abs(x);

  if (abs_x >= 0.0 && abs_x <= 1.0)
  {
    return 1.0 - (abs_x / 2.0) - pow(abs_x, 2.0) + (pow(abs_x, 3.0) / 2.0);
  }
  else if (abs_x >= 1.0 && abs_x <= 2.0)
  {
    return 1.0 - ((11.0 * abs_x) / 6.0) + pow(abs_x, 2.0) - (pow(abs_x, 3.0) / 6.0);
  }

  return 0;
}

double w_6(double x)
{
  auto abs_x = abs(x);

  if (abs_x >= 0.0 && abs_x <= 1.0)
  {
    return 1.0 - (abs_x / 3.0) - ((5.0 * pow(abs_x, 2.0)) / 4.0) + ((5.0 * pow(abs_x, 3.0)) / 12.0) + (pow(abs_x, 4.0) / 4.0) - (pow(abs_x, 5.0) / 12.0);
  }
  else if (abs_x >= 1.0 && abs_x <= 2.0)
  {
    return 1.0 - ((13.0 * abs_x) / 12.0) - ((5.0 * pow(abs_x, 2.0)) / 8.0) + ((25.0 * pow(abs_x, 3.0)) / 24.0) - ((3.0 * pow(abs_x, 4.0)) / 8.0) + (pow(abs_x, 5.0) / 24.0);
  }
  else if (abs_x >= 2.0 && abs_x <= 3.0)
  {
    return 1.0 - ((137.0 * abs_x) / 60.0) + ((15.0 * pow(abs_x, 2.0)) / 8.0) - ((17.0 * pow(abs_x, 3.0)) / 24.0) + (pow(abs_x, 4.0) / 8.0) - (pow(abs_x, 5.0) / 120.0);
  }

  return 0;
}

double w_helper(array<double, DIM> z, const double& h_g, double (*w)(double))
{
  double product = 1;
  for (int i=0; i<DIM; ++i)
  {
    product *= (1./h_g)*w(z[i]/h_g);
  }
  return product;
}

w_ptr get_w(const int& choice)
{
  switch(choice)
  {
    case 2:
      return w_2;
    case 3:
      return w_3;
    case 4:
      return w_4;
    case 6:
      return w_6;
    default:
      throw new invalid_argument("Not a valid choice for spline interpolant!");
  }
}

double w(array<double, DIM> z, double h_g, int choice)
{
  double result = 0;
  switch (choice)
  {
    case 2:
      return w_helper(z, h_g, w_2);
    case 3:
      return w_helper(z, h_g, w_3);
    case 4:
      return w_helper(z, h_g, w_4);
    case 6:
      return w_helper(z, h_g, w_6);
    default:
      throw new invalid_argument("Not a valid choice for spline interpolant!");
  }
}