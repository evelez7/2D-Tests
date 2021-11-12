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

  return (a*d) - (b*c);
}

array<array<double, DIM>, DIM> multiply_matrices(const array<array<double, DIM>, DIM>& first, const array<array<double, DIM>, DIM>& second)
{
  array<array<double, DIM>, DIM> product;
  for (int i=0; i<DIM;++i)
  {
    for (int j=0; j<DIM;++j)
    {
      product[i][j] = 0.;
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

bool is_identity(const matrix &m, const double& margin)
{
  if (abs(m[0][0] - 1.) >= margin)
    return false;
  if (abs(m[0][1] - 0.) >= margin)
    return false;
  if (abs(m[1][0] - 0.) >= margin)
    return false;
  if (abs(m[1][1] - 1.) >= margin)
    return false;
  return true;
}

// check if matrices are equal according to a margin of error between respective indices
bool matrices_are_equal(const matrix &m_1, const matrix &m_2, const double &margin)
{
  for (int i=0; i<DIM; ++i)
    for (int j=0; j<DIM; ++j)
    {
      if (abs(m_1[i][j] - m_2[i][j]) <= margin)
        continue;
      else
        return false;
    }
  return true;
}

// detailed check on whether or not R and Q are valid matrices
// Q will be an orthogonal rotation matrix
void verify_R_Q(const matrix &R, const matrix& Q, const matrix& deformation_matrix)
{
  cout << boolalpha;
  cout << endl << "NEW VERIFICATION" << endl;
  bool R_inverse_valid = is_identity(multiply_matrices(R, get_inverse(R)), 0.001);
  bool R_transpose_valid =matrices_are_equal(get_transpose(R), R);
  auto R_R = multiply_matrices(R, R);
  bool R_R_valid = matrices_are_equal(R_R, multiply_matrices(get_transpose(deformation_matrix), deformation_matrix), 0.1);
  cout << "verify for R" << endl;
  cout << "R^T = R : " <<  R_transpose_valid << endl;
  cout << "R * R^-1 = I : " << R_inverse_valid << endl;
  cout << "R * R = A^T * A: " << R_R_valid << endl;
  cout << endl;

  cout << "verify for Q" << endl;
  bool Q_inverse_transpose = matrices_are_equal(get_inverse(Q), get_transpose(Q), 0.001);
  bool Q_identity = is_identity(multiply_matrices(Q, get_transpose(Q)), 0.001);
  bool Q_R_mult = matrices_are_equal(multiply_matrices(Q, R), deformation_matrix, 0.001);
  cout << "Q^-1 = Q^T : " <<  Q_inverse_transpose << endl;
  cout << "Q * Q^T = I : " <<  Q_identity << endl;
  cout << "Q * R = A : " <<  Q_R_mult << endl;
  cout << endl;

  if (!R_inverse_valid || !R_transpose_valid || !Q_inverse_transpose || Q_identity || Q_R_mult)
  {

  cout << "A" << endl;
  print_matrix(deformation_matrix);
  cout << "Q * R = A" << endl;
  print_matrix(multiply_matrices(Q, R));
  cout << "R_R" << endl;
  print_matrix(R_R);
  cout << "R" << endl;
  print_matrix(R);
  cout << "R^T" << endl;
  print_matrix(get_transpose(R));
  cout << "R^-1" << endl;
  print_matrix(get_inverse(R));
  cout << "Q" << endl;
  print_matrix(Q);
  cout << "Q^-1" << endl;
  print_matrix(get_inverse(Q));
  cout << "Q^T" << endl;
  print_matrix(get_transpose(Q));
  }

  bool q_det;
  if (get_determinant(Q) == 1.)
    q_det = true;
  else
    q_det = false;
  cout << "det(Q) = 1 : " << q_det << " it is " << get_determinant(Q) << endl;
  cout << endl;
}

bool matrices_are_equal(const array<array<double, DIM>, DIM>& matrix_1, const array<array<double, DIM>, DIM>& matrix_2)
{
  for (int i=0; i<DIM; ++i)
    for (int j=0; j<DIM; ++j)
      if (matrix_1[i][j] != matrix_2[i][j]) return false;

  return true;
}

matrix subtract_matrices(const matrix& m_1, const matrix& m_2)
{
  matrix to_return;
  for (int i=0; i<DIM; ++i)
  {
    for (int j=0; j<DIM; ++j)
    {
      to_return[i][j] = m_1[i][j] - m_2[i][j];
    }
  }
  return to_return;
}

matrix get_identity_matrix()
{
  matrix identity_matrix;
  identity_matrix[0][0] = 1.;
  identity_matrix[0][1] = 0.;
  identity_matrix[1][0] = 0.;
  identity_matrix[1][1] = 1.;
  return identity_matrix;
}