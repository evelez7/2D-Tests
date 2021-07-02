#ifndef VELOCITIES_H
#define VELOCITIES_H

#include "Proto_Point.H"
#include <array>

using namespace std;
using namespace Proto;

double find_magnitude(const Point&);

double find_magnitude(const array<double, DIM>&);

double find_velocity(const double&, const int&);

double find_velocity_derivative(const double&, const double&);

double f(const double&);

#endif