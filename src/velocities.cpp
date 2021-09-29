#include <cmath>
#include <exception>
#include "velocities.h"

/**
 * @brief calculate angular velocity formula
 *
 * @param radius the distance of the particle from the center
 * @param p commandline input
 * @return double velocity of the particle
 */
double find_velocity(const double &r, const int &p, const double &r0)
{
  if (r < r0)
    return (1. - pow((1. - pow(r / r0, 2)), p + 1)) / (2. * static_cast<double>(p + 1) * pow(r / r0, 2));
  else if (r > r0)
    return 1. / (2. * static_cast<double>(p + 1) * pow(r / r0, 2));
  else if (r0 == 0.)
    return 0.;
  throw new runtime_error("Reached end of find_velocity without falling into an if clause");
}

// TODO
double find_velocity_derivative(const double &r, const int &p, const double &r0)
{
  double p_d = static_cast<double>(p);
  // derivative with respect to r
  if (r < r0)
  {
    double numerator = (pow(r,2)*pow((pow(r0,2)-pow(r,2)), p)*pow(r0,-2*(p+1))+pow(r,2)*p*pow((pow(r0,2)-pow(r,2)), p)*pow(r0, -2*(p+1))+pow((-(pow(r,2)/pow(r0,2))+1),p+1))*pow(r0,2);
    double denominator = (p+1)*pow(r,3);
    return numerator / denominator;
  }
  else if (r > r0)
    return -((pow(r0,2))/(pow(r, 3)*(p+1)));
  else if (r0==0.)
    return 0.;
  throw new runtime_error("Reached end of find_velocity_derivative without falling into an if clause");
}

double find_magnitude(const array<double, DIM> &alpha)
{
  return sqrt(pow(alpha[0], 2) + pow(alpha[1], 2));
}

double find_magnitude(const double &alpha_d)
{
  return sqrt(pow(alpha_d, 2.));
}

double f(const double &r, const double &r0)
{
  if (r <= r0)
    return 1.;
  else if (r > r0)
    return 0.;
  throw new runtime_error("Reached end of exact function without falling into an if clause");
}