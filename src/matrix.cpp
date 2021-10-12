#include "matrix.h"
#include "Proto_Point.H"

array<array<double, DIM>, DIM> get_transpose(const array<array<double, DIM>, DIM>& matrix) {
  array<array<double, DIM>, DIM> transpose;

  for (int i=0; i < DIM; ++i)
  {
    for (int j=0; j < DIM; ++j)
    {
      transpose[i][j] = matrix[j][i];
    }
  }
  return transpose;
}

double get_determinant(const array<array<double, DIM>, DIM>& matrix)
{
  double a = matrix.at(0).at(0);
  double b = matrix.at(0).at(1);
  double c = matrix.at(1).at(0);
  double d = matrix.at(1).at(1);

  return (a*d) + (-1. * (b*c));
}

array<array<double, DIM>, DIM> multiply_matrices(const array<array<double, DIM>, DIM>& first, const array<array<double, DIM>, DIM>& second)
{
  array<array<double, DIM>, DIM> product;
  for (int i=0; i<DIM;++i)
  {
    for (int j=0; j<DIM;++j)
    {
      product[i][j] = 0;
      for (int k=0; k<DIM;++k)
      {
        product[i][j] += first[i][k] * second[k][j];
      }
    }
  }
  return product;
}

array<double, DIM> multiply_matrix_by_vector(const array<array<double, DIM>, DIM>& matrix, const array<double, DIM>& vec)
{
  array<double, DIM> new_vector;
  for (int i=0;i<DIM;++i) {
    new_vector[i] = 0;
    for (int j=0; j < matrix[i].size(); ++j) {
      new_vector[i] +=  matrix[i][j] * vec[j];
    }
  }
  return new_vector;
}

// only works for
array<array<double, DIM>, DIM> get_inverse(const array<array<double, DIM>, DIM>& matrix)
{
  array<array<double, DIM>, DIM> matrix_inverse;
  double determinant = 1./((matrix[0][0]*matrix[1][1]) - (matrix[0][1]*matrix[1][0]));
  matrix_inverse[0][0] = matrix[1][1] * determinant;
  matrix_inverse[0][1] = -matrix[0][1] * determinant;
  matrix_inverse[1][0] = -matrix[1][0] * determinant;
  matrix_inverse[1][1] = matrix[0][0] * determinant;
  return matrix_inverse;
}

void print_matrix(const array<array<double, DIM>, DIM>& matrix)
{
  for (int i=0; i<DIM; ++i)
  {
    for (int j=0; j<DIM; ++j)
    {
      cout << matrix[i][j] << " ";
    }
    cout << endl;
  }
  cout << endl;
}

void print_vector(const array<double, DIM> &vec)
{
  cout << "(";
  for (int i=0; i<DIM; ++i)
    cout << vec[i] << ", ";
  cout << ")" << endl;
}

array<double, DIM> multiply_vector_by_scalar(const array<double, DIM>& vec, const double& scalar)
{
  array<double, DIM> result;
  for (int i=0; i<DIM; ++i)
    result.at(i) = vec.at(i) * scalar;
  return result;
}

matrix multiply_matrix_by_scalar(matrix to_multiply, const double& scalar)
{
  for (int i=0; i<DIM; ++i)
    for (int j=0; j<DIM; ++j)
      to_multiply[i][j] = to_multiply[i][j] * scalar;
  return to_multiply;
}

/**
 * @brief Get the unit vector in the d direction
 *
 * @param d
 * @return array<double, DIM>
 */
array<double, DIM> get_unit_vector(const int& d)
{
  array<double, DIM> unit_vector;
  for (int i=0; i<d; ++i)
  {
    if (i == d-1)
    {
      unit_vector[i] = 1.;
      continue;
    }
    unit_vector[i] = 0.;
  }
  return unit_vector;
}

matrix get_zero_matrix()
{
  matrix zero;
  zero[0][0] = 0.;
  zero[0][1] = 0.;
  zero[1][0] = 0.;
  zero[1][1] = 0.;
  return zero;
}