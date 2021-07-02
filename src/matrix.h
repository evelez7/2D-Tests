#ifndef MATRIX_H
#define MATRIX_H

#include <array>
#include <cmath>
#include <iostream>

using namespace std;

typedef array<array<double, DIM>, DIM> matrix;

array<array<double, DIM>, DIM> get_transpose(const array<array<double, DIM>, DIM>&);

double get_determinant(const array<array<double, DIM>, DIM>&);

array<array<double, DIM>, DIM> multiply_matrices(const array<array<double, DIM>, DIM>&, const array<array<double, DIM>, DIM>&);

array<double, DIM> multiply_matrix_by_vector(const array<array<double, DIM>, DIM>&, const array<double, DIM>&);

array<array<double, DIM>, DIM> get_inverse(const array<array<double, DIM>, DIM>&);

void print_matrix(const array<array<double, DIM>, DIM>&);

array<double, 3> get_characteristic_polynomial(const array<array<double, DIM>, DIM>&);

array<double, DIM> multiply_vector_by_scalar(const array<double, DIM>&, const double&);

matrix multiply_matrix_by_scalar(matrix, const double&);

array<double, DIM> get_unit_vector(const int&);

#endif