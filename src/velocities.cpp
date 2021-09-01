#include <cmath>
#include <exception>
#include "velocities.h"

double r0;

void initialize_r0(const double& new_r0)
{
  r0 = new_r0;
}

/**
 * @brief calculate angular velocity formula
 *
 * @param radius the distance of the particle from the center
 * @param p commandline input
 * @return double velocity of the particle
 */
double find_velocity(const double& r, const int& p)
{

  if (r <= r0)
    return static_cast<double>((1. - pow(1. - static_cast<double>(pow(r/r0, 2.)), p+1)))/(2. * (p+1) * static_cast<double>(pow(r/r0, 2.)));
  else if (r >= r0)
    return 1./(2. * static_cast<double>(p+1.) * static_cast<double>(pow(r/r0, 2.)));

}

// TODO
double find_velocity_derivative(const double& r, const int& p)
{
  double p_d = static_cast<double>(p);
  // derivative with respect to r
  // if (r <= r0)
  // {
  //   double numerator = ((pow(r,2.) * p_d) + pow(r0, 2.)) *pow((pow(r0, 2.) - pow(r, 2.)), p_d);
  //   double denominator = pow(r, 3.)*pow(r0, 2.*p_d)*(p_d+1);
  //   return numerator/denominator;
  // } else if (r >= r0)
  // {
  //   return -pow(r0, 2.)/(pow(r, 3.) * (p +1.));
  // }

  // with respect to p
  if ( r <= r0)
  {
    double lhs_num = pow(1. - (pow(r, 2.)/pow(r0, 2.)), p_d+1.) * log(1.- (pow(r, 2.)/pow(r0,2.))) * (p_d +1.);
    double rhs_num = pow(1. - (pow(r, 2.)/pow(r0, 2.)), p_d+1.);
    double num = (lhs_num - rhs_num) * pow(r0, 2.);

    double denom = 2. * pow(r, 2.) * pow(p+1, 2.);

    return num/denom;
  } else if (r >= r0)
  {
    return -pow(r0, 2.)/(2. * pow(r, 2.) * pow(p+1., 2.));
  }

  // return 1.;
}

double find_magnitude(const array<double, DIM>& alpha)
{
  double magnitude = sqrt(pow(alpha[0], 2.) + pow(alpha[1], 2.));
  return magnitude;
}

double find_magnitude(const double& alpha_d)
{
  return sqrt(pow(alpha_d, 2.));
}

double f(const double& r)
{
  if (r <= r0) return 1.;
  else if (r >= r0) return 0.;
  return 0.;
}