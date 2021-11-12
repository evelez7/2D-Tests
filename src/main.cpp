#include <cmath>
#include <array>
#include <vector>
#include <iostream>
#include <string>
#include <memory>
#include "matrix.h"
#include "w.h"
#include "eigen_util.cpp"
#include "velocities.h"
#include "writers.h"
#include "interpolate.h"
#include "particles.h"
#include "Proto_Box.H"
#include "Proto_BoxData.H"
#include "Proto_Point.H"
#include "Proto_WriteBoxData.H"

static FILE *fp = NULL;
static int useBinary = 0;
static int numInColumn = 0;

using namespace std;
using namespace Proto;

array<double, DIM> lagrangian_map(const array<double, DIM> &, const double &, const int &, const double &);

matrix rotation(const array<double, DIM> &, const double &, const int &, const double &);

matrix partial_derivative_rotation(const array<double, DIM> &, const array<double, DIM> &, const double &, const int &, const int &, const double &);

matrix compute_deformation_matrix(const Particle &, const Particle &, const double &, const double &, const double &);
vector<matrix> compute_deformation_matrices(const vector<Particle> &, const vector<Particle> &, const double &, const double &, const double &);

array<matrix, 2> compute_R_Q(const matrix &, const array<double, DIM>&);

void print_array(const array<double, DIM> &to_print)
{
  cout << "{ ";
  for (int i = 0; i < DIM; ++i)
    cout << to_print[i] << " ";
  cout << "}" << endl;
}

// initialize particles spaced by h_p
vector<Particle> initialize_particles(const double &hp, const int &np)
{
  vector<Particle> particles;
  for (int i = static_cast<int>(-np / 2); i < static_cast<int>(np / 2); ++i)
  {
    for (int j = static_cast<int>(-np / 2); j < static_cast<int>(np / 2); ++j)
    {
      particles.emplace_back(static_cast<double>(i) * hp, static_cast<double>(j) * hp, pow(hp, 2.), 0., 0., 0., 0.);
    }
  }
  return particles;
}

vector<Particle> move_particles(const vector<Particle> particles, const double &time, const int &p, const double& r0)
{
  vector<Particle> rotated_particles;
  for (int i = 0; i < particles.size(); ++i)
  {
    array<double, DIM> part_coords{particles[i].x, particles[i].y};
    auto projection = lagrangian_map(part_coords, time, p, r0);
    rotated_particles.emplace_back(projection[0], projection[1], particles[i].strength, particles[i].velocity, 0., 0., 0.);
  }
  return rotated_particles;
}

int main(int argc, char **argv)
{
  int p, log_n, np, ng;

  cout << "number of grid points (log_2):" << endl;
  cin >> log_n;
  ng = static_cast<int>(pow(2, log_n));

  cout << "Enter p:" << endl;
  cin >> p;

  double time;
  cout << "Enter time:" << endl;
  cin >> time;
  // time *= M_PI;

  cout << endl
       << "Spline choices are: " << endl;
  for (int i = 1; i <= get_interpolant_size(); ++i)
  {
    cout << i << ".) " << get_interpolant_name(i) << endl;
  }

  int spline_choice;
  cout << "Enter spline interpolant:" << endl;
  cin >> spline_choice;
  int ppercell;
  cout << "enter number of particles per cell:" << endl;
  cin >> ppercell;
  double r0;
  cout << "enter r0:" << endl;
  cin >> r0;
  np = ng * ppercell;
  Point bottom_left(-ng, -ng);
  Point top_right(ng, ng);
  Box grid_box(bottom_left, top_right);
  BoxData<double> grid(grid_box);
  grid.setVal(0.);

  double hg = 1. / static_cast<double>(ng); // interparticle spacing
  double hp = hg / static_cast<double>(ppercell);

  vector<Particle> particles = initialize_particles(hp, np);
  vector<Particle> rotated_particles = move_particles(particles, time, p, r0);
  for (int i = 0; i < rotated_particles.size(); ++i)
  {
    array<double, DIM> x_k {rotated_particles[i].x, rotated_particles[i].y};
    interpolate(grid, rotated_particles[i].strength, x_k, hg, hp, spline_choice);

    array<array<double, DIM>, DIM> deformation_matrix = compute_deformation_matrix(particles[i], rotated_particles[i], time, p, r0);
    auto A_t_A = multiply_matrices(get_transpose(deformation_matrix), deformation_matrix); // the symmetric and positive definite matrix
    //equation 30-32
    auto eigenvalues = find_trace_eigenvalues(A_t_A);
    auto grad_det = get_determinant(A_t_A);
    double eigen_product = eigenvalues[0] * eigenvalues[1];

    rotated_particles[i].eigen_1 = eigenvalues[0];
    rotated_particles[i].eigen_2 = eigenvalues[1];
    rotated_particles[i].eigen_product = eigen_product;
    auto r_q = compute_R_Q(deformation_matrix, eigenvalues);
    verify_R_Q(r_q[0], r_q[1], deformation_matrix);
  }
  string filename = "grid";
  cout << "Writing grid to files" << endl;
  WriteData(grid, 0, hg, filename, filename);
  cout << "Writing particles to file" << endl;
  PWrite(rotated_particles);
  cout << "Finished" << endl;

  return 0;
}

// return a deformation matrix from individual particles
matrix compute_deformation_matrix(const Particle &particle, const Particle &rotated_particle, const double &time, const double &p, const double &r0)
{
  matrix deformation_matrix;
  deformation_matrix[0][0] = 0.;
  deformation_matrix[0][1] = 0.;
  deformation_matrix[1][0] = 0.;
  deformation_matrix[1][1] = 0.;
  // fill deformation matrix
  array<double, DIM> alpha{particle.x, particle.y};
  array<double, DIM> x_k{rotated_particle.x, rotated_particle.y};
  for (int j = 0; j < DIM; ++j)
  {
    auto lhs = multiply_matrix_by_vector(rotation(alpha, time, p, r0), get_unit_vector(j + 1));
    auto rhs = multiply_matrix_by_vector(partial_derivative_rotation(x_k, alpha, time, p, j, r0), alpha);
    if (j == 0)
    {
      deformation_matrix[0][0] = lhs[0] + rhs[0];
      deformation_matrix[1][0] = lhs[1] + rhs[1];
    }
    else if (j == 1)
    {
      deformation_matrix[0][1] = lhs[0] + rhs[0];
      deformation_matrix[1][1] = lhs[1] + rhs[1];
    }
  }
  return deformation_matrix;
}

// returns a vector of deformation matrices in order of the original particle vector
vector<matrix> compute_deformation_matrices(const vector<Particle> &particles, const vector<Particle> &rotated_particles, const double &time, const double &p, const double &r0)
{
  vector<matrix> deformation_matrices;
  for (int i = 0; i < particles.size(); ++i)
  {
    array<double, DIM> alpha{particles[i].x, particles[i].y};
    array<double, DIM> x_k{rotated_particles[i].x, rotated_particles[i].y};

    array<array<double, DIM>, DIM> deformation_matrix;
    deformation_matrix[0][0] = 0.;
    deformation_matrix[0][1] = 0.;
    deformation_matrix[1][0] = 0.;
    deformation_matrix[1][1] = 0.;
    // fill deformation matrix
    for (int j = 0; j < DIM; ++j)
    {
      auto lhs = multiply_matrix_by_vector(rotation(alpha, time, p, r0), get_unit_vector(j + 1));
      auto rhs = multiply_matrix_by_vector(partial_derivative_rotation(x_k, alpha, time, p, j, r0), alpha);

      if (j == 0)
      {
        deformation_matrix[0][0] = lhs[0] + rhs[0];
        deformation_matrix[1][0] = lhs[1] + rhs[1];
      }
      else if (j == 1)
      {
        deformation_matrix[0][1] = lhs[0] + rhs[0];
        deformation_matrix[1][1] = lhs[1] + rhs[1];
      }
    }
    deformation_matrices.push_back(deformation_matrix);
  }
  return deformation_matrices;
}

// returns an array where arr[0] = R, arr[1] = Q
array<matrix, 2> compute_R_Q(const matrix &deformation_matrix, const array<double, DIM>& eigenvalues)
{
  auto A_t_A = multiply_matrices(get_transpose(deformation_matrix), deformation_matrix); // the symmetric and positive definite matrix
  // equation 30-32
  auto E = find_trace_eigenvectors(A_t_A, eigenvalues);
  auto E_transposed = get_transpose(E);

  array<array<double, DIM>, DIM> diag;
  diag[0][0] = eigenvalues[0];
  diag[0][1] = 0;
  diag[1][1] = eigenvalues[1];
  diag[1][0] = 0;

  auto R_lhs = multiply_matrices(E, diag);
  auto R = multiply_matrices(R_lhs, E_transposed); // symm matrix
  // auto R = multiply_matrices(E_transposed, R_lhs);
  auto Q = multiply_matrices(deformation_matrix, get_inverse(R)); // rotation
  array<matrix, 2> matrices;
  matrices[0] = R;
  matrices[1] = Q;

  return matrices;
}

array<double, DIM> lagrangian_map(const array<double, DIM> &alpha, const double &time, const int &p, const double &r0)
{
  auto rotation_matrix = rotation(alpha, time, p, r0);
  auto projection = multiply_matrix_by_vector(rotation_matrix, alpha);
  return projection;
}

matrix rotation(const array<double, DIM> &alpha, const double &time, const int &p, const double &r0)
{
  double velocity = 0.;
  if (find_magnitude(alpha) != 0.)
  {
    velocity = find_velocity(find_magnitude(alpha), p, r0);
  }
  matrix rotation_matrix;
  rotation_matrix[0][0] = cos(velocity * time);
  rotation_matrix[0][1] = sin(velocity * time);
  rotation_matrix[1][0] = -sin(velocity * time);
  rotation_matrix[1][1] = cos(velocity * time);

  return rotation_matrix;
}

/**
 * @brief The partial derivative of the rotation matrix with respect to the Lagrangian coordinate
 *
 * @param alpha
 * @param time
 * @param p
 * @return matrix
 */
matrix partial_derivative_rotation(const array<double, DIM> &alpha, const array<double, DIM> &original_alpha, const double &time, const int &p, const int &d, const double &r0)
{
  // edge case for the very center particle
  if (find_magnitude(original_alpha) == 0.)
  {
    return get_zero_matrix();
  }
  double velocity = find_velocity(find_magnitude(original_alpha), p, r0);
  // partial derivative with respect to alpha_d
  double velocity_derivative = find_velocity_derivative(find_magnitude(original_ialpha), p, r0);
  double velocity_deriv_alpha = velocity_derivative * (original_alpha[d] / find_magnitude(alpha));
  matrix rotation_matrix;
  rotation_matrix[0][0] = -sin(velocity_deriv_alpha * time) * time * velocity_deriv_alpha;
  rotation_matrix[0][1] = cos(velocity_deriv_alpha * time) * time * velocity_deriv_alpha;
  rotation_matrix[1][0] = -cos(velocity_deriv_alpha * time) * time * velocity_deriv_alpha;
  rotation_matrix[1][1] = -sin(velocity_deriv_alpha * time) * time * velocity_deriv_alpha;

  return rotation_matrix;
}
