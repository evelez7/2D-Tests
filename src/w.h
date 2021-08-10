#ifndef W_H
#define W_H

#include <array>
#include <string>

using namespace std;

typedef double (*w_ptr)(double);

double w_2(double);

double w_3(double);

double w_4(double);

double w_6(double);

double L2_1(double);

double L2_2(double);

double L4_2(double);

double L4_4(double);

/**
 * @brief choose which spline interpolant to use
 *
 * @return double
 */
double w(array<double, DIM>, double, int);

w_ptr get_w(const int&);

string get_interpolant_name(const int&);

int get_interpolant_size();

#endif