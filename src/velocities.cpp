#include <cmath>
#include <exception>
#include "velocities.h"

double r0=1;

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
    return (1 - pow(1 - pow(r/r0, 2.), p+1))/(2 * (p+1) * pow(r/r0, 2.));
  else if (r >= r0)
    return 1/(2 * (p+1) * pow(r/r0, 2.));
    // return 1;
}

double find_velocity_derivative(const double& radius, const int& p)
{

}

double find_magnitude(const Point& alpha)
{
  // double alpha_1 = abs(static_cast<double>(alpha[0]) - center[0]);
  // double alpha_2 = abs(static_cast<double>(alpha[1]) - center[1]);
  double magnitude = sqrt(pow(alpha[0], 2.) + pow(alpha[1], 2.));
  // cout << "mag: " << magnitude << " x: " << alpha_1 << " y: " << alpha_2 << endl;
  return magnitude;
}

double find_magnitude(const array<double, DIM>& alpha)
{
  double magnitude = sqrt(pow(alpha[0], 2.) + pow(alpha[1], 2.));
  // cout << "mag: " << magnitude << endl;
  return magnitude;
}

double f(const double& r)
{
  if (r <= r0) return 1.;
  else if (r >= r0) return 0.;
}