#include <cmath>
#include <array>
#include <vector>
#include <iostream>
#include <string>
#include <memory>
#include "matrix.h"
#include "w.h"
#include "velocities.h"
#include "interpolate.h"
#include "Proto_Box.H"
#include "Proto_BoxData.H"
#include "Proto_Point.H"
#include "Proto_WriteBoxData.H"
#include "Proto_VisitWriter.H"
static FILE *fp = NULL;
static int useBinary = 0;
static int numInColumn = 0;

using namespace std;
using namespace Proto;

struct Particle
{
  double x;
  double y;
  double strength;
  double velocity;
  Particle(double x, double y, double strength, double velocity)
        : x(x)
        , y(y)
        , strength(strength)
        , velocity(velocity)
    {}
};

array<double, DIM> lagrangian_map(const array<double, DIM>&, const double&, const int&);

/**
 * @brief Get the rotation for a specific coordinate at time t
 *
 * @return matrix
 */
matrix rotation(const Point &, const double &, const int &);

matrix rotation(const array<double, DIM>&, const double &, const int &);

array<double, DIM> partial_derivative_x(const Point &, const double &, const int &, const int &);

array<double, DIM> init_point(const Point& j, const double &h_p, const double &time, const int &p)
{
  array<double, DIM> j_scaled {{ j[0] * h_p, j[1] * h_p }};
  return multiply_matrix_by_vector(rotation(j_scaled, time, p), j_scaled);
}

inline void write_point_mesh(const char *filename, int npts, double *pts,
                             int nvars, int *vardim, const char *const *varnames,
                             double **vars)
{
  FILE *fp = vtk_open_file(filename);
  int numInColumn = 0;

  int i;
  char str[128];
  int *centering = NULL;

  write_header(fp);

  write_string("DATASET UNSTRUCTURED_GRID\n", fp);
  sprintf(str, "POINTS %d double\n", npts);
  write_string(str, fp);
  for (i = 0; i < 3 * npts; i++)
  {
    write_double(pts[i], fp, numInColumn);
  }

  new_section(fp, numInColumn);
  sprintf(str, "CELLS %d %d\n", npts, 2 * npts);
  write_string(str, fp);
  for (i = 0; i < npts; i++)
  {
    write_int(1, fp, numInColumn);
    write_int(i, fp, numInColumn);
    end_line(fp, numInColumn);
  }
  centering = (int *)malloc(nvars * sizeof(int));
  for (i = 0; i < nvars; i++)
    centering[i] = 1;
  write_variables(nvars, vardim, centering, varnames, vars, npts, npts, fp, numInColumn);
  free(centering);

  vtk_close_file(fp, numInColumn);
}

void PWrite(const char *a_filename, const vector<Particle> particles, int fileCount)
{
  vector<vector<double>> vars(4);
  unsigned int size = particles.size();
  vector<double> x(3 * size);
  vars[0] = vector<double>(size);
  vars[1] = vector<double>(size);
  vars[2] = vector<double>(size);
  vars[3] = vector<double>(size);
  for (unsigned int i = 0; i < size; ++i)
  {
    auto current_part = particles.at(i);
    vars[0][i] = current_part.x;
    vars[1][i] = current_part.y;
    vars[2][i] = current_part.strength;
    vars[3][i] = current_part.velocity;
    x[i * 3] = current_part.x;
    x[i * 3 + 1] = current_part.y;
    x[i * 3 + 2] = 0.;
    // cout << "x: " << current_part.x << " y: " << current_part.y << endl;
  }
  double *varPtr[4];
  varPtr[0] = &vars[0][0];
  varPtr[1] = &vars[1][0];
  varPtr[2] = &vars[2][0];
  varPtr[3] = &vars[3][0];
  int vardim[4] = {1, 1, 1, 1};
  const char *const varnames[] = {"x_1", "x_2", "strength", "velocity"};
  write_point_mesh("parts", size, &(x[0]), 4, vardim, varnames, varPtr);
}

void PWrite(vector<Particle> particles)
{
  static int fileCount = 0;
  static char nameBuffer[10];
  sprintf(nameBuffer, "PART.%d", fileCount);
  PWrite(nameBuffer, particles, fileCount);
}

void print_array(const array<double, DIM>& to_print)
{
  cout << "{ ";
  for (int i=0; i<DIM; ++i)
    cout << to_print[i] << " ";
  cout << "}" << endl;
}

// initialize particles spaced by h_p
vector<Particle> initialize_particles(const double& hp, const int& np)
{
  array<double, DIM> new_coords;
  vector<Particle> particles;
  for (int i=-np/2; i<np/2; ++i)
  {
    new_coords[0] = i*hp;
    for (int j=-np/2; j<np/2; ++j)
    {
      new_coords[1] = j*hp;
      Particle new_part { new_coords[0], new_coords[1], pow(hp, 2.), 0. };
      particles.push_back(new_part);
    }
  }
  return particles;
}

vector<Particle> move_particles(vector<Particle> particles, const double& time, const int& p)
{
  vector<Particle> rotated_particles;
  for (int i=0; i<particles.size(); ++i)
  {
    array<double, DIM> part_coords { particles[i].x, particles[i].y };
    if (part_coords[0] == 0 && part_coords[1] == 0)
    {
      rotated_particles.emplace_back(0., 0., particles[i].strength, particles[i].velocity);
      continue;
    }
    auto projection = lagrangian_map(part_coords, time, p);
    rotated_particles.emplace_back(projection[0], projection[1], particles[i].strength, particles[i].velocity);
  }
  return rotated_particles;
}

int main(int argc, char **argv)
{
  int p, log_n, np,ng;

  cout << "number of grid points (log_2): ";
  cin >> log_n;
  ng = static_cast<int>(pow(2, log_n));

  cout << "Enter p: ";
  cin >> p;
  // array<double, DIM> center = { static_cast<double>(range)/2., static_cast<double>(range)/2. };

  double time;
  cout << "Enter time: ";
  cin >> time;
  time *= M_PI;
  int spline_choice;
  cout << "Enter spline interpolant: ";
  cin >> spline_choice;
  int ppercell;
  cout << "enter number of particles per cell" << endl;
  cin >> ppercell;
  // Box grid_box(Point::Zeros(), Point::Unit()*np); // coordinates from (0,0) to (p, p)
  np=ng*ppercell;
  Point bottom_left(-ng, -ng);
  Point top_right(ng, ng);
  Box grid_box(bottom_left, top_right);
  BoxData<double> grid(grid_box);

  // double hp = 1.; // interparticle spacing
  // double hg = 1.;
  double hg = 1./static_cast<double>(ng); // interparticle spacing
  double hp = hg/ppercell;
  cout << "Particles per cell: " << hg*hg/hp/hp << endl;

  auto particles = initialize_particles(hp, np);
  // vector<Particle> particles = initialize_particles(grid.box(), hp, np);
  vector<Particle> rotated_particles = move_particles(particles, time, p);

  for (int i=0; i<rotated_particles.size(); ++i)
  {
    array<double, DIM> x_k {rotated_particles[i].x, rotated_particles[i].y};
    interpolate(grid, rotated_particles[i].strength, x_k, hg, hp, spline_choice);
  }
  string filename = "grid";
  WriteData(grid, 0, hg, filename, filename);
  PWrite(rotated_particles);

  return 0;
}

array<double, DIM> lagrangian_map(const array<double, DIM>& alpha, const double& time, const int& p)
{
  return multiply_matrix_by_vector(rotation(alpha, time, p), alpha);
}

matrix rotation(const array<double, DIM> &alpha, const double&time, const int &p)
{
  double velocity = find_velocity(find_magnitude(alpha), p);
  // cout << "velocity: " << velocity << endl;
  matrix rotation_matrix;
  rotation_matrix[0][0] = cos(velocity * static_cast<double>(time) );
  rotation_matrix[0][1] = sin(velocity * static_cast<double>(time) );
  rotation_matrix[1][0] = -sin(velocity * static_cast<double>(time) );
  rotation_matrix[1][1] = cos(velocity * static_cast<double>(time) );

  return rotation_matrix;
}

matrix rotation(const Point &alpha, const double&time, const int &p)
{
  array<double, DIM> alpha_array { static_cast<double>(alpha[0]), static_cast<double>(alpha[1]) };
  double velocity = find_velocity(find_magnitude(alpha_array), p);
  matrix rotation_matrix;
  rotation_matrix[0][0] = cos(velocity * static_cast<double>(time) );
  rotation_matrix[0][1] = sin(velocity * static_cast<double>(time) );
  rotation_matrix[1][0] = -sin(velocity * static_cast<double>(time) );
  rotation_matrix[1][1] = cos(velocity * static_cast<double>(time) );

  return rotation_matrix;
}

// matrix partial_derivative_velocity(const double&

/**
 * @brief The partial derivative of the rotation matrix with respect to the Lagrangian coordinate
 *
 * @param alpha
 * @param time
 * @param p
 * @return matrix
 */
matrix partial_derivative_rotation(const Point &alpha, const double &time, const int &p)
{
  double velocity = find_velocity(find_magnitude(alpha), p);
  matrix rotation_matrix;
  rotation_matrix[0][0] = -sin(velocity * time) * time;
  rotation_matrix[0][1] = cos(velocity * time) * time;
  rotation_matrix[1][0] = -cos(velocity * time) * time;
  rotation_matrix[1][1] = -sin(velocity * time) * time;
}

/**
 * @brief
 *
 * @param alpha
 * @param time
 * @param d
 * @param p
 * @return array<double, DIM>
 */
array<double, DIM> partial_derivative_x(const Point &alpha, const double &time, const int &d, const int &p)
{
  auto left_hand_side = multiply_matrix_by_vector(rotation(alpha, time, p), get_unit_vector(d));

  array<double, DIM> alpha_array{{static_cast<double>(alpha[0]), static_cast<double>(alpha[1])}};
  auto right_hand_side = multiply_matrix_by_vector(partial_derivative_rotation(alpha, time, p), alpha_array);

  array<double, DIM> to_return{{left_hand_side[0] + left_hand_side[1], right_hand_side[0] + right_hand_side[1]}};
  return to_return;
}
