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
static int printDone = 0;
using namespace std;
using namespace Proto;

array<double, DIM> lagrangian_map(const array<double, DIM> &, const double &, const int &, const double &);

matrix rotation(const array<double, DIM> &, const double &, const int &, const double &);

matrix partial_derivative_rotation(const array<double, DIM> &, const array<double, DIM> &, const double &, const int &, const int &, const double &);

matrix compute_deformation_matrix(const Particle &, const Particle &, const double &, const double &, const double &);
vector<matrix> compute_deformation_matrices(const vector<Particle> &, const vector<Particle> &, const double &, const double &, const double &);

array<matrix, 2> compute_R_Q(const matrix &);
vector<Particle> work(int, int, int, int, double, double);
void test_strengths_plots(vector<Particle>, string);
void test_minmax(vector<Particle>, string, double);

int main(int argc, char** argv)
{
  int p_minmax = 0;
  int p_strengths = 8;
  vector<int> ppercell_vec = {1, 2, 4};
  double r0=0.5;
  double log2 = 6;
  for (int i=1; i <= get_interpolant_size(); ++i)
    for (auto ppercell : ppercell_vec)
    {
      string param_string = to_string(ppercell) + "_" + get_interpolant_name(i) + "_r" + to_string(r0) + "_log" + to_string(log2) + ".curve";
      // i will be spline choice
      // last param is time
      // ppercell, interpolant, p, log2(n), r0, time
      auto gridded_particles_strengths = work(ppercell, i, p_strengths, log2, r0, 0.5);
      test_strengths_plots(gridded_particles_strengths, param_string);

      // for the min max graphs, do separately
      string minmax_filename = "./tests/min_max/min_max_pc" + param_string;
      fstream minmax_curve_file;

      minmax_curve_file.open(minmax_filename, ios::app | ios::out);
      minmax_curve_file << "# " << minmax_filename << endl;
      minmax_curve_file.close();
      for (double t=0.05; t<=0.5; t+=0.05)
      {
        auto gridded_particles_minmax = work(ppercell, i, p_minmax, log2, r0, t);
        test_minmax(gridded_particles_minmax, param_string, t);
        cout << "TEST" << endl;
      }
    }

  return 0;
}

void print_array(const array<double, DIM> &to_print)
{
  cout << "{ ";
  for (int i = 0; i < DIM; ++i)
    cout << to_print[i] << " ";
  cout << "}" << endl;
}
void print_matrix(array<array<double, DIM>, DIM> a_mat,string a_str)
{
  if (printDone == 0)
    {
      cout << "matrix " << a_str << endl;
      cout << a_mat[0][0] << " , " << a_mat[0][1] << endl;
      cout << a_mat[1][0] << " , " << a_mat[1][1] << endl << endl;
    }
}
bool check_orthogonal(matrix Q)
{
  double eps=1.e-7;
  auto idprod = multiply_matrices(Q,get_transpose(Q));
  bool ret = ((fabs(idprod[0][0] - 1.) < eps) && (fabs(idprod[1][1] - 1.) < eps) && (fabs(idprod[0][1]) < eps) && (fabs(idprod[1][0]) < eps));
  return ret;
}
bool check_eigenvectors(matrix& E,matrix& A, double lambda0,double lambda1)
{
  double eps = 1.e-7;
  // Eigenvectors are the columns of E.
  auto AE = multiply_matrices(A,E);
  bool ret = fabs(AE[0][0] - lambda0*E[0][0]) < max(fabs(E[1][0]),fabs(E[0][0]))*eps;
  //cout << "ret = " << ret << endl;
  if (!ret)
    {
    cout << AE[0][0] - 1 << " , " << lambda0*E[0][0] - 1 << " ," << fabs(E[1][0]) << " , " << fabs(E[0][0]) << endl;
    cout << "eigenvalues =" << lambda0 << " , " << lambda1 << endl;
    }
  ret = ret && fabs(AE[1][0] - lambda0*E[1][0]) < max(fabs(E[1][0]),fabs(E[0][0]))*eps;
  //cout << "ret = " << ret << endl;
  ret = ret && fabs(AE[0][1] - lambda1*E[0][1]) < max(fabs(E[0][1]),fabs(E[1][1]))*eps;
  //cout << "ret = " << ret << endl;
  ret = ret && fabs(AE[1][1] - lambda1*E[1][1]) < max(fabs(E[0][1]),fabs(E[1][1]))*eps;
  //cout << "ret = " << ret << endl;
  return ret;
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

vector<Particle> work(int ppercell, int spline_choice, int p, int log_n, double r0, double time)
{
  int np, ng;
  // p=8;log_n=6;time= .5*M_PI;spline_choice = 6;r0=.5;ppercell=2;
  time *= M_PI;

  ng = static_cast<int>(pow(2, log_n));
  np = ng * ppercell;
  Point bottom_left(-ng, -ng);
  Point top_right(ng, ng);
  Box grid_box(bottom_left, top_right);
  BoxData<double> gridded_angle(grid_box);
  BoxData<double> gridded_eigen(grid_box);
  BoxData<double> gridded_strength(grid_box);
  gridded_angle.setVal(0.);
  gridded_eigen.setVal(0.);
  gridded_strength.setVal(0.);

  double hg = 1. / static_cast<double>(ng); // interparticle spacing
  double hp = hg / static_cast<double>(ppercell);

  vector<Particle> particles = initialize_particles(hp, np);
  vector<Particle> rotated_particles = move_particles(particles, time, p, r0);

  for (int i = 0; i < rotated_particles.size(); ++i)
  {
    array<double, DIM> x_k{rotated_particles[i].x, rotated_particles[i].y};

    array<array<double, DIM>, DIM> deformation_matrix =
      compute_deformation_matrix(particles[i], rotated_particles[i], time, p, r0);
    auto A_t_A =
      multiply_matrices(get_transpose(deformation_matrix), deformation_matrix);
    // the symmetric and positive definite matrix
    //equation 30-32
    auto eigenvalues = get_sym_eigenvalues(A_t_A);
    auto grad_det = get_determinant(A_t_A);
    auto rq = compute_R_Q(deformation_matrix);
    double detQ = get_determinant(rq[1]);
    double eigen_product = eigenvalues[0] * eigenvalues[1];
    rotated_particles[i].angle = acos(rq[1][0][0]);
    rotated_particles[i].eigen_1 = eigenvalues[0];
    rotated_particles[i].eigen_2 = eigenvalues[1];
    interpolate(gridded_angle, rotated_particles[i].angle*hp*hp, x_k, hg, hp, spline_choice);
    interpolate(gridded_eigen, rotated_particles[i].eigen_1*hp*hp, x_k, hg, hp, spline_choice);
    interpolate(gridded_strength, rotated_particles[i].strength, x_k, hg, hp, spline_choice);
  }

  {
    string filename = "gridded_angle";
    cout << "Writing gridded data to files" << endl;
    WriteData(gridded_angle, 0, hg, filename, filename);
  }
  {
    string filename = "gridded_eigen";
    WriteData(gridded_eigen, 0, hg, filename, filename);
  }
  {
    string filename = "gridded_strength";
    WriteData(gridded_strength, 0, hg, filename, filename);
  }
  cout << "Writing particles to file" << endl;
  PWrite(rotated_particles);
  std::ofstream ofs ("outparticles1D.txt", std::ofstream::out);

  cout << "Generate new particles corresponding to interpolated values, write to file" << endl;
  vector<Particle> gridded_particles;
  vector<Particle> old_gridded_particles;
  for (auto bit = grid_box.begin(); !bit.done(); ++bit)
    {
      Point pos = *bit;
      array<double,DIM> rpos;
      rpos[0] = pos[0]*hg;
      rpos[1] = pos[1]*hg;
      double rmag =max(abs(rpos[0]),abs(rpos[1]));
      double gridStrength = gridded_strength(pos);
      Particle pt(rpos[0],rpos[1],gridStrength,0,0,0,0);
      if (rmag < .3) gridded_particles.push_back(pt);
    }

  old_gridded_particles = move_particles(gridded_particles,-time,p,r0);
  for (int i = 0; i < gridded_particles.size();i++)
    {
      array<array<double, DIM>, DIM> deformation_matrix = compute_deformation_matrix(old_gridded_particles[i], gridded_particles[i], time, p, r0);

      auto A_t_A = multiply_matrices(get_transpose(deformation_matrix), deformation_matrix);
      auto eigenvalues = get_sym_eigenvalues(A_t_A);
      auto grad_det = get_determinant(A_t_A);
      double eigen_product = eigenvalues[0] * eigenvalues[1];
      auto rq = compute_R_Q(deformation_matrix);
      double detQ = get_determinant(rq[1]);
      auto alpha0 = old_gridded_particles[i].x;
      auto alpha1 = old_gridded_particles[i].y;
      array<double,DIM> alpha = {alpha0,alpha1};
      auto amag = sqrt(pow(alpha0,2) + pow(alpha1,2));
      array<double,DIM> ahat = {alpha0/amag,alpha1/amag};
      array<double,DIM> aperp = {ahat[1],-ahat[0]};
      array<array<double,DIM>,DIM> EalphaT = {ahat,aperp};
      auto Ealpha = get_transpose(EalphaT);
      gridded_particles[i].eigen_1 = eigenvalues[0];
      auto R = rotation(alpha, time, p, r0);
      auto Qinvariant =
      multiply_matrices(Ealpha,
                        multiply_matrices(R,multiply_matrices(rq[1],EalphaT)));
      gridded_particles[i].eigen_2 = eigenvalues[1];
      gridded_particles[i].angle = acos(rq[1][0][0]);
      if (!check_orthogonal(rq[1]))
        {
          cout << "orthogonal matrix check fail: rq[1] =" << endl;
          print_matrix(rq[1]);
          abort();
        }
      ofs << i << " " << eigenvalues[0] << " " << eigen_product << " " << rq[1][0][0] << " "<< gridded_particles[i].strength - 1.0 << endl;
      //ofs << grad_det << " " << gridded_particles[i].strength - 1.0 << endl;
    }
  string fn = "gridded_particles";
  string ofn = "old_gridded_particles";
  PWrite(fn.c_str(),gridded_particles,0);
  PWrite(ofn.c_str(),old_gridded_particles,0);
  ofs.close();

  cout << "Finished" << endl;
  return gridded_particles;
}

void test_strengths_plots(vector<Particle> gridded_particles, string param_string)
{
    vector<double> strength_axis, eigen_1_axis, angle_axis;
    for (auto particle : gridded_particles)
    {
      strength_axis.push_back(particle.strength-1.);
      eigen_1_axis.push_back(particle.eigen_1);
      angle_axis.push_back(particle.angle);
    }
    string strength_eigen1_filename ="./tests/strengths/eigen/strength_eigen1_pc";
    write_curve_file(strength_axis, eigen_1_axis, strength_eigen1_filename + param_string);
    string strength_angle_filename ="./tests/strengths/angle/strength_angle_pc";
    write_curve_file(strength_axis, angle_axis, strength_angle_filename + param_string);
    cout << "Done writing strength vs y tests..." << endl;
}

void test_minmax(vector<Particle> gridded_particles, string filename_params, double time)
{
    double min = numeric_limits<double>::max();
    double max = numeric_limits<double>::min();
    for (auto particle : gridded_particles)
    {
      if (particle.strength - 1. > max)
        max = particle.strength;
      if (particle.strength - 1. < min)
        min = particle.strength;
    }
    vector<double> x,y;
    x = {time, time};
    y = {min, max};
    cout << min << " " << max << endl;
    string filename = "./tests/min_max/min_max_pc" + filename_params;
    write_curve_file_append(x, y, filename);
}

// return a deformation matrix from individual particles
matrix compute_deformation_matrix(const Particle &particle, const Particle &rotated_particle, const double &time, const double &p, const double &r0)
{
  matrix deformation_matrix;
  deformation_matrix[0][0] = 1.;
  deformation_matrix[0][1] = 0.;
  deformation_matrix[1][0] = 0.;
  deformation_matrix[1][1] = 1.;
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
array<matrix, 2> compute_R_Q(const matrix &deformation_matrix)
{
  auto A_t_A = multiply_matrices(get_transpose(deformation_matrix), deformation_matrix);
  // the symmetric and positive definite matrix
  // equation 30-32
  auto eigenvalues = get_sym_eigenvalues(A_t_A);
  auto grad_det = get_determinant(A_t_A);
  double eigen_product = eigenvalues[0] * eigenvalues[1];
  if (fabs(sqrt(grad_det)-eigen_product) > 1.e-12)
    {
      cout << "grad_det = " << grad_det << " , eigen_product = " << eigen_product << endl;
    }
  auto E_transposed = find_sym_eigenvectors(A_t_A, eigenvalues);
  // conveniently, this is E transposed
  auto E = get_transpose(E_transposed); // do not confuse eigenvectors with E

  auto E_inverse = get_inverse(E);
  array<array<double, DIM>, DIM> diag;
  diag[0][0] = eigenvalues[0];
  diag[0][1] = 0;
  diag[1][1] = eigenvalues[1];
  diag[1][0] = 0;
  auto chk_eigen = check_eigenvectors(E,A_t_A, pow(eigenvalues[0],2), pow(eigenvalues[1],2));
  if (!chk_eigen)
    {
      cout << "not eigenvectors!"<< endl;
      auto AE = multiply_matrices(A_t_A,E);
      cout << "AE = " << endl;
      print_matrix(AE);
      cout << "E = " << endl;
      print_matrix(E);
      cout << "eigenvalues - 1 = " <<  pow(eigenvalues[0],2) - 1 << " , " <<  pow(eigenvalues[1],2) - 1 << endl;
      abort();
    }
  auto R = multiply_matrices(multiply_matrices(E, diag),E_inverse); // symm matrix
  auto R_inverse = get_inverse(R);
  auto Q = multiply_matrices(deformation_matrix, R_inverse); // rotation
  array<matrix, 2> matrices;
  matrices[0] = R;
  matrices[1] = Q;
  // check correctness
  if (!check_orthogonal(Q))
    {
      cout << "Q not orthogonal, det(Q) = "<< get_determinant(Q) << endl;
      //print_matrix(Q);
      //print_matrix(R);
      print_matrix(multiply_matrices(R,R));
      print_matrix(A_t_A);
      auto AE = multiply_matrices(A_t_A,E);
      print_matrix(AE);
      print_matrix(multiply_matrices(diag,E));
      abort();
    }
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
  double magnitude = find_magnitude(alpha);
  double velocity = 0.;
  if (magnitude != 0.)
    velocity = find_velocity(magnitude, p, r0);
  matrix rotation_matrix;
  rotation_matrix[0][0] = cos(velocity * time);
  rotation_matrix[0][1] = sin(velocity * time);
  rotation_matrix[1][0] = -sin(velocity * time);
  rotation_matrix[1][1] = cos(velocity * time);
  //print_matrix(rotation_matrix);

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
  double velocity_derivative = find_velocity_derivative(find_magnitude(original_alpha), p, r0);
  double velocity_deriv_alpha = velocity_derivative * (original_alpha[d] / find_magnitude(original_alpha));
  matrix rotation_matrix;
  rotation_matrix[0][0] = -sin(velocity * time) * time * velocity_deriv_alpha;
  rotation_matrix[0][1] = cos(velocity * time) * time * velocity_deriv_alpha;
  rotation_matrix[1][0] = -cos(velocity * time) * time * velocity_deriv_alpha;
  rotation_matrix[1][1] = -sin(velocity * time) * time * velocity_deriv_alpha;

  return rotation_matrix;
}