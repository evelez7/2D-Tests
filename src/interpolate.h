#ifndef INTERPOLATE_H
#define INTERPOLATE_H

#include <array>
#include <vector>
#include "Proto_BoxData.H"

using namespace std;

// array<std::array<double, DIM>, DIM> interpolate(const Proto::BoxData<double> [DIM][DIM], const Particle&, const double&);

// array<double, DIM> interpolate(const array<double, DIM>&, const array<double, DIM>&, const array<double, DIM>, const double&, const int&);

double interpolate(const double&, const array<double, DIM>, const double&, const int&);


#endif