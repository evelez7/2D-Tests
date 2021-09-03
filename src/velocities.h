#ifndef VELOCITIES_H
#define VELOCITIES_H

#include "Proto_Point.H"
#include <array>

using namespace std;
using namespace Proto;

double find_magnitude(const array<double, DIM>&);

double find_velocity(const double&, const int&, const double&);

double find_velocity_derivative(const double&, const int&, const double&);

double f(const double&, const double&);

#endif