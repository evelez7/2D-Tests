#ifndef PARTICLES_H
#define PARTICLES_H

struct Particle
{
  double x;
  double y;
  double strength;
  double velocity;
  double eigen_1;
  double eigen_2;
  double eigen_product;
  Particle(double x, double y, double strength, double velocity, double eigen_1, double eigen_2, double eigen_product)
      : x(x), y(y), strength(strength), velocity(velocity), eigen_1(eigen_1), eigen_2(eigen_2), eigen_product(eigen_product)
  {
  }
};

#endif