#ifndef RESERVE_H_INCLUDED
#define RESERVE_H_INCLUDED
#include <iostream>
#include <math.h>
using namespace std;

class reserve
{
private:
  float x, y;
  int dimension;
  friend float operator+(const reserve &r1, const reserve &r2);

public:
  void setx(float);
  void sety(float);
  float getx() const;
  float gety() const;
  int getDim() const;
  float operator[](const int index) const;
  float addReserve(reserve const &r2) const;
};

#endif
