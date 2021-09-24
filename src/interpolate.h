#ifndef INTERPOLATE_H
#define INTERPOLATE_H

#include <array>
#include <vector>
#include "Proto_BoxData.H"

using namespace std;
using namespace Proto;

void interpolate(BoxData<double>&, const double&, const array<double, DIM>, const double&, const double&, const int&);

void interpolate_test(const double&, const array<double, DIM>, const double&, const double&, const int&);

#endif
