#include "matrix.h"
#include <iostream>
#include <cmath>
using namespace std;

// arr[0] = a, arr[1] = b, arr[2] = c, x=-b +- sqrt(b^2 -4ac) all over 2a
array<double, DIM> get_roots(const array<double, 3>& characteristic_polynomial)
{
  double a = characteristic_polynomial[0];
  double b = characteristic_polynomial[1];
  double c = characteristic_polynomial[2];
  double discriminant = pow(b, 2.) - (4. * a * c);
  array<double, DIM> roots;
  if (discriminant > 0)
  {
    roots[0] = (-b + sqrt(discriminant))/(2*a);
    roots[1] = (-b - sqrt(discriminant))/(2*a);
  } else if (discriminant == 0 || discriminant < 1e-10)
  {
    roots[0] = -b/(2.*a);
    roots[1]=roots[0];
  }
  // if (discriminant < 1e-10)
  // {
  //   cout << "root 1: " << roots[0] << endl;
  //   cout << "root 2: " << roots[1] << endl;
  // }
  return roots;
}

array<double, 3> get_characteristic_polynomial(const array<array<double, DIM>, DIM>& matrix)
{
  auto right_det = (matrix[0][1] * matrix[1][0]); // the rhs of determinant, b*c
  auto constant = matrix[0][0] * matrix[1][1];
  auto lambda_coeff_1 = -matrix[0][0];
  auto lambda_coeff_2 = -matrix[1][1];

  array<double, 3> characteristic_polynomial;
  characteristic_polynomial[0] = 1;
  characteristic_polynomial[1] = lambda_coeff_1 + lambda_coeff_2;
  characteristic_polynomial[2] = constant - right_det;
  return characteristic_polynomial;
}

bool is_zero_matrix(const array<array<double, DIM>, DIM>& matrix)
{
  int threshold = 4;
  int count = 0;
  for (auto row : matrix)
    for (auto elem : row)
      if (elem == 0)
        return true;
  return false;
}

void verify_eigenvalues(const array<array<double, DIM>, DIM>& matrix, const array<double, DIM>& eigen_values)
{
  cout << "DETERMINANTS: ";
  for (auto lambda : eigen_values)
  {
    array<array<double, DIM>, DIM> to_check;
    to_check[0][0] = matrix[0][0] - lambda;
    to_check[0][1] = matrix[0][1];
    to_check[1][0] = matrix[1][0];
    to_check[1][1] = matrix[1][1] - lambda;

    double determinant = (to_check[0][0] * to_check[1][1]) - (to_check[0][1]*to_check[1][0]);
    cout << determinant << " ";
  }
  cout << endl;
}

// prints matrix for the eigenvectors, and the eigenvalues and eigenvectors
void verify_eigenvectors_verbose(const array<array<double, DIM>, DIM>& matrix, const array<array<double, DIM>, DIM>& eigenvectors, const array<double, DIM>& eigenvalues)
{
  cout << "VERIFYING EIGENVECTORS FOR MATRIX" << endl;
  print_matrix(matrix);
  for (int i=0; i<DIM; ++i)
  {
      cout << "EIGENVALUE " << i << ": " << eigenvalues[i] <<  endl;
      auto rhs = multiply_vector_by_scalar(eigenvectors[i], eigenvalues[i]);
      auto lhs = multiply_matrix_by_vector(matrix, eigenvectors[i]);
      cout << "LHS" << endl << "A * " << eigenvalues[i] << endl << "[";
      for (auto elem : lhs)
        cout << elem << " ";
      cout << "]" << endl;
      cout << "RHS" << endl <<  "eigenvector: [";
      for (auto elem : eigenvectors[i])
        cout << elem << " ";
      cout << "] * " << eigenvalues[i] << endl << "[";
      for (auto elem : rhs)
        cout << elem << " ";
      cout << "]" << endl << endl;
  }
  cout << endl;
}

// push back errors of the components of the eigenvectors to later see the max and min errors between the components
void collect_eigenvector_errors(const array<array<double, DIM>, DIM>& matrix, const array<array<double, DIM>, DIM>& eigenvectors, const array<double, DIM>& eigenvalues, shared_ptr<vector<double>>& collection)
{
  for (int i=0; i<DIM; ++i)
  {
      auto rhs = multiply_vector_by_scalar(eigenvectors.at(i), eigenvalues.at(i));
      auto lhs = multiply_matrix_by_vector(matrix, eigenvectors.at(i));
      for (int j=0; j<DIM; ++j)
      {
        collection->push_back(lhs.at(j) - rhs.at(j));
      }
  }
}

// visually verify eigenvector accuracy by printing the error between the components
void verify_eigenvectors(const array<array<double, DIM>, DIM>& matrix, const array<array<double, DIM>, DIM>& eigenvectors, const array<double, DIM>& eigenvalues)
{
  cout << "verifying eigenvectors of matrix" << endl;
  print_matrix(matrix);
  for (int i=0; i<DIM; ++i)
  {
      auto rhs = multiply_vector_by_scalar(eigenvectors[i], eigenvalues[i]);
      auto lhs = multiply_matrix_by_vector(matrix, eigenvectors[i]);
      for (int j=0; j<DIM; ++j)
      {
        double error = lhs.at(j) - rhs.at(j);
        cout << "error: " << error << " ";
      }
      cout << endl;
  }
  cout << endl;
}

array<array<double, DIM>, DIM> find_eigenvectors(const array<array<double, DIM>, DIM> &A, const array<double, DIM> &eigenvalues)
{
  // take the special case that we have a 2x2 matrix
  array<array<double, DIM>, DIM> eigenvectors;

  double a = A[0][0];
  double b = A[0][1];
  double c = A[1][0];
  double d = A[1][1];

  if ((c == 0 || c < 1e-5) && (b == 0 || b < 1e-5))
  {
    eigenvectors[0][0] = 1;
    eigenvectors[0][1] = 0;
    eigenvectors[1][0] = 0;
    eigenvectors[1][1] = 1;
  } else if (c != 0 || c > 1e-5)
  {
    eigenvectors[0][0] = eigenvalues[0] - d;
    eigenvectors[0][1] = c;
    eigenvectors[1][0] = eigenvalues[1] - d;
    eigenvectors[1][1] = c;
  } else if (b != 0 || b > 1e-5)
  {
    eigenvectors[0][0] = b;
    eigenvectors[0][1] = eigenvalues[0] - a;
    eigenvectors[1][0] = b;
    eigenvectors[1][1] = eigenvalues[1] - a;
  }

  return eigenvectors;
}
array<array<double, DIM>, DIM> find_sym_eigenvectors(const array<array<double, DIM>, DIM> &A, const array<double, DIM> &eigenvalues)
{
  // take the special case that we have a symmetric 2x2 matrix
  array<array<double, DIM>, DIM> eigenvectors;  

  double eps = 1.e-14;
  double b = A[1][0];
  if (abs(b - A[0][1]) > eps) abort();
  double a = A[0][0] - pow(eigenvalues[0],2);
  double d = A[1][1] - pow(eigenvalues[0],2);
  int l0,l1;
  if (abs(a) > abs(d)) l0 = 0; else l0 = 1;
  l1 = (l0 + 1) % DIM;
  a = A[l0][l0] - pow(eigenvalues[0],2);
  // The eigenvector for lambda is orthogonal to the rows of ATilde = A - lambda*I.
  // For the 2x2 case, we can compute the eigenvector as ATilde_{l,.}^\perp.
  if ((abs(b) > abs(a))&&(abs(a) > eps))
    {
      eigenvectors[0][l0] = 1.;
      eigenvectors[0][l1] = -a/b;
    }
  else if ((abs(a) > abs(b))&&(abs(b) > eps))
    {
      eigenvectors[0][l1] = -1.;
      eigenvectors[0][l0] = b/a;
    }
  else
    {
      eigenvectors[0][l0] = 1.;
      eigenvectors[0][l1] = 0.;
    }
  // For a symmetric matrix, the eigenvectors are orthogonal.
  eigenvectors[1][l0] = -eigenvectors[0][l1];
  eigenvectors[1][l1] = eigenvectors[0][l0];
  return eigenvectors;
}
// equation 30-32
// returns the square root of the eigenvalues of the symmetric matrix R of the polar decomposition of the gradient
array<double, DIM> get_sym_eigenvalues(const array<array<double, DIM>, DIM>& matrix)
{
  // equation 31
  auto poly = get_characteristic_polynomial(matrix);
  auto roots = get_roots(poly);

  // visual verification where the determinants should be close to 0
  // verify_eigen_values(A_t_A, roots);

  array<double, DIM> eigen_diag;
  // equation 32
  for (int i=0; i<DIM; ++i)
    eigen_diag[i] = sqrt(roots[i]);

  return eigen_diag;
}

array<double, DIM> get_eigenvalues(const array<array<double, DIM>, DIM>& matrix)
{

  // equation 31
  auto poly = get_characteristic_polynomial(matrix);
  auto roots = get_roots(poly);

  // visual verification where the determinants should be close to 0
  // verify_eigen_values(A_t_A, roots);

  array<double, DIM> eigen_diag;
  // equation 32
  for (int i=0; i<DIM; ++i)
    eigen_diag[i] = roots[i];

  return eigen_diag;
}

array<double, DIM> get_eigenvalues_trace(const array<array<double, DIM>, DIM>& matrix)
{
  double trace = matrix[0][0] + matrix[1][1];
  double determinant = (matrix[0][0] * matrix[1][1]) - (matrix[0][1] * matrix[1][0]);
  array<double, DIM> eigenvalues;
  eigenvalues[0] = (trace/2) + sqrt( (pow(trace, 2.)/4) - determinant);
  eigenvalues[1] = (trace/2) - sqrt( (pow(trace, 2.)/4) - determinant);

  return eigenvalues;
}

double get_magnitude(const array<double, DIM> &vec)
{
  double sum=0;
  for (auto elem : vec)
    sum += pow(elem, 2);
  return sqrt(sum);
}

bool matrices_are_equal(const array<array<double, DIM>, DIM>& matrix_1, const array<array<double, DIM>, DIM>& matrix_2)
{
  for (int i=0; i<DIM; ++i)
    for (int j=0; j<DIM; ++j)
      if (matrix_1[i][j] != matrix_2[i][j]) return false;

  return true;
}
