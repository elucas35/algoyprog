#ifndef LEAF_H_INCLUDED
#define LEAF_H_INCLUDED
#include <vector>
#include "reserve.h"
class leaf{
private:
  vector< reserve > coordleaf;
  float minx;
  float maxx;
  float miny;
  float maxy;

public:
  void addToLeaf(float* coord);
  vector<reserve> getCoordLeaf();
  void setMinx(float minx);
  void setMiny(float miny);
  void setMaxx(float maxx);
  void setMaxy(float maxy);
  float getMinx();
  float getMiny();
  float getMaxx();
  float getMaxy();
};
#endif
