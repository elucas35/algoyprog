#include "leaf.h"
void leaf::addToLeaf(float* coord){
 reserve r;
 r.setx(coord[0]);
 r.sety(coord[1]);
 coordleaf.push_back(r);
}
vector<reserve> leaf::getCoordLeaf(){
 return coordleaf;
}
void leaf::setMinx(float min){
  minx = min;
}
void leaf::setMiny(float min){
  miny = min; 
}
void leaf::setMaxx(float max){
  maxx = max;
}
void leaf::setMaxy(float max){
  maxy = max;
}
float leaf::getMinx(){
  return minx;
}
float leaf::getMiny(){
  return miny;
}
float leaf::getMaxx(){
  return maxx;
}
float leaf::getMaxy(){
  return maxy;
}

