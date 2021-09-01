#ifndef VELOCITIES_H
#define VELOCITIES_H

#include "Proto_Point.H"
#include <array>

using namespace std;
using namespace Proto;

void initialize_r0(const double&);

double find_magnitude(const Point&);

double find_magnitude(const array<double, DIM>&);

double find_velocity(const double&, const int&);

double find_velocity_derivative(const double&, const int&);

double f(const double&);

#endif