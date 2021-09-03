#ifndef WRITERS_H
#define WRITERS_H

#include <vector>
#include "particles.h"

using namespace std;

inline void write_point_mesh(const char *, int, double *, int , int *, const char *const *, double **);

void PWrite(const char *, const vector<Particle>, int);

void PWrite(vector<Particle>);

#endif