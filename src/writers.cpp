#include "Proto_WriteBoxData.H"
#include "Proto_VisitWriter.H"
#include "writers.h"

void PWrite(vector<Particle> particles)
{
  static int fileCount = 0;
  static char nameBuffer[10];
  sprintf(nameBuffer, "PART.%d", fileCount);
  PWrite(nameBuffer, particles, fileCount);
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
  vector<vector<double>> vars(7);
  unsigned int size = particles.size();
  vector<double> x(2 * size);
  vars[0] = vector<double>(size);
  vars[1] = vector<double>(size);
  vars[2] = vector<double>(size);
  vars[3] = vector<double>(size);
  vars[4] = vector<double>(size);
  vars[5] = vector<double>(size);
  vars[6] = vector<double>(size);
  for (unsigned int i = 0; i < size; ++i)
  {
    auto current_part = particles.at(i);
    vars[0][i] = current_part.x;
    vars[1][i] = current_part.y;
    vars[2][i] = current_part.strength;
    vars[3][i] = current_part.velocity;
    vars[4][i] = current_part.eigen_1;
    vars[5][i] = current_part.eigen_2;
    vars[6][i] = current_part.eigen_product;

    x[i * 3] = current_part.x;
    x[i * 3 + 1] = current_part.y;
    x[i * 3 + 2] = 0.;
    // cout << "x: " << current_part.x << " y: " << current_part.y << endl;
  }
  double *varPtr[7];
  varPtr[0] = &vars[0][0];
  varPtr[1] = &vars[1][0];
  varPtr[2] = &vars[2][0];
  varPtr[3] = &vars[3][0];
  varPtr[4] = &vars[4][0];
  varPtr[5] = &vars[5][0];
  varPtr[6] = &vars[6][0];
  int vardim[7] = {1, 1, 1, 1, 1, 1, 1};
  const char *const varnames[] = {"x_1", "x_2", "strength", "velocity", "eigen_1", "eigen_2", "eigen_product"};
  write_point_mesh("parts", size, &(x[0]), 7, vardim, varnames, varPtr);
}