#include <cmath>
#include <iostream>
#include <stdexcept>
#include "w.h"

using namespace std;

array<string, 8> interpolant_names{
    "w_2", "w_3", "w_4", "w_6", "L2_1", "L2_2", "L4_2", "L4_4"};

double w_2(double x)
{
  double abs_x = abs(x);
  if (abs_x >= 0. && abs_x <= 1.)
  {
    return 1. - abs_x;
  }
  return 0.;
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

  return 0.;
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

  return 0.;
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

  return 0.;
}

double L2_1(double x)
{
  if (0 <= abs(x) && abs(x) < 1.)
    return 1. - ((5. / 2.) * (pow(abs(x), 2))) + ((3. / 2.) * (pow(abs(x), 3)));
  else if (1 <= abs(x) && abs(x) < 2.)
    return 2. - (4. * (abs(x))) + ((5. / 2.) * (pow(abs(x), 2))) - ((pow(abs(x), 3.)) / 2.);
  else if (2 <= abs(x))
    return 0.;
}

double L2_2(double x)
{
  if (0. <= abs(x) && abs(x) < 1.)
    return 1. - (pow(abs(x), 2)) - ((9. / 2.) * (pow(abs(x), 3))) + ((15. / 2.) * (pow(abs(x), 4))) - (3. * (pow(abs(x), 5)));
  else if (1. <= abs(x) && abs(x) < 2.)
    return (-4. + (18. * abs(x)) - (29. * (pow(abs(x), 2))) + ((43. / 2.) * (pow(abs(x), 3))) - ((15. / 2.) * (pow(abs(x), 4))) + (pow(abs(x), 5)));
  else if (2 <= abs(x))
    return 0.;
}

double L4_2(double x)
{
  if (0 <= abs(x) && abs(x) < 1.)
    return (1. - ((5. / 4.) * (pow(abs(x), 2))) - ((35. / 12.) * (pow(abs(x), 3))) + ((21. / 4.) * (pow(abs(x), 4))) - ((25. / 12.) * (pow(abs(x), 5))));
  else if (1. <= abs(x) && abs(x) < 2.)
    return (-4. + ((75. / 4.) * (abs(x))) - ((245. / 8.) * (pow(abs(x), 2))) + ((545. / 24.) * (pow(abs(x), 3))) - ((63. / 8.) * (pow(abs(x), 4))) + ((25. / 24.) * (pow(abs(x), 5))));
  else if (2. <= abs(x) && abs(x) < 3.)
    return (18. - ((153. / 4.) * (abs(x))) + ((255. / 8.) * (pow(abs(x), 2))) - ((313. / 24.) * (pow(abs(x), 3))) + ((21. / 8.) * (pow(abs(x), 4))) - ((5. / 24.) * (pow(abs(x), 5))));
  else if (3 <= abs(x))
    return 0.;
}

double L4_4(double x)
{
  if (0 <= abs(x) && abs(x) < 1)
    return (1. - ((5. / 4.) * pow(abs(x), 2)) + ((1. / 4.) * pow(abs(x), 4)) - ((100. / 3.) * pow(abs(x), 5)) + ((455. / 4.) * pow(abs(x), 6)) - ((295. / 2.) * pow(abs(x), 7)) + ((345. / 4.) * pow(abs(x), 8)) - ((115. / 6.) * pow(abs(x), 9)));
  else if (1 <= abs(x) && abs(x) < 2.)
    return (-199. + ((5485. / 4.) * (abs(x))) - ((32975. / 8.) * pow(abs(x), 2)) + ((28425. / 4.) * pow(abs(x), 3)) - ((61953. / 8.) * pow(abs(x), 4)) + ((33175. / 6.) * pow(abs(x), 5)) - ((20685. / 8.) * pow(abs(x), 6)) + ((3055. / 4.) * pow(abs(x), 7)) - ((1035. / 8.) * pow(abs(x), 8)) + ((115. / 12.) * pow(abs(x), 9)));
  else if (2. <= abs(x) && abs(x) < 3)
    return (5913. - ((89235. / 4.) * (abs(x))) + ((297585. / 8.) * pow(abs(x), 2)) - ((143895. / 4.) * pow(abs(x), 3)) + ((177871. / 8.) * pow(abs(x), 4)) - ((54641. / 6.) * pow(abs(x), 5)) + ((19775. / 8.) * pow(abs(x), 6)) - ((1715. / 4.) * pow(abs(x), 7)) + ((345. / 8.) * pow(abs(x), 8)) - ((23. / 12.) * pow(abs(x), 9)));
  else if (abs(x) >= 3)
    return 0;
}

double w_helper(array<double, DIM> z, const double &h_g, w_ptr w)
{
  double product = 1.;
  for (int i = 0; i < DIM; ++i)
  {
    product *= (1. / h_g) * w(z[i] / h_g);
  }
  return product;
}

w_ptr get_w(const int &choice)
{
  switch (choice)
  {
  case 1:
    return w_2;
  case 2:
    return w_3;
  case 3:
    return w_4;
  case 4:
    return w_6;
  case 5:
    return L2_1;
  case 6:
    return L2_2;
  case 7:
    return L4_2;
  case 8:
    return L4_4;
  default:
    throw new invalid_argument("Not a valid choice for spline interpolant!");
  }
}

string get_interpolant_name(const int &choice)
{
  return interpolant_names[choice - 1];
}

int get_interpolant_size()
{
  return interpolant_names.size();
}

double w(array<double, DIM> z, double h_g, int choice)
{
  switch (choice)
  {
  case 1:
    return w_helper(z, h_g, w_2);
  case 2:
    return w_helper(z, h_g, w_3);
  case 3:
    return w_helper(z, h_g, w_4);
  case 4:
    return w_helper(z, h_g, w_6);
  case 5:
    return w_helper(z, h_g, L2_1);
  case 6:
    return w_helper(z, h_g, L2_2);
  case 7:
    return w_helper(z, h_g, L4_2);
  case 8:
    return w_helper(z, h_g, L4_4);
  default:
    throw new invalid_argument("Not a valid choice for spline interpolant!");
  }
}