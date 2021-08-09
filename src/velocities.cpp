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

  // if (r <= r0)
    return static_cast<double>((1. - pow(1. - static_cast<double>(pow(r/r0, 2.)), p+1)))/(2. * (p+1) * static_cast<double>(pow(r/r0, 2.)));
  // else if (r >= r0)
  //   return 1./(2. * static_cast<double>(p+1.) * static_cast<double>(pow(r/r0, 2.)));

}

// TODO
double find_velocity_derivative(const double& radius, const int& p)
{

}

double find_magnitude(const Point& alpha)
{
  double magnitude = sqrt(pow(alpha[0], 2.) + pow(alpha[1], 2.));
  return magnitude;
}

double find_magnitude(const array<double, DIM>& alpha)
{
  double magnitude = sqrt(pow(alpha[0], 2.) + pow(alpha[1], 2.));
  return magnitude;
}

double f(const double& r)
{
  if (r <= r0) return 1.;
  else if (r >= r0) return 0.;
}