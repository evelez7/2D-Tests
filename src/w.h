#ifndef W_H
#define W_H

#include <array>

using namespace std;

typedef double (*w_ptr)(double);

double w_2(double);

double w_3(double);

double w_4(double);

double w_6(double);

/**
 * @brief choose which spline interpolant to use
 *
 * @return double
 */
double w(array<double, DIM>, double, int);

w_ptr get_w(const int&);

#endif