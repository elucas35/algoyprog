
#include "reserve.h"
void reserve::setx(float coordx){
  x = coordx;
}
void reserve::sety(float coordy){
  y = coordy;
}
float reserve::getx() const{
  return x;
}
float reserve::gety() const{
  return y;
}
int reserve::getDim() const{
  cout << "2" << endl;
  return 2;
}
float reserve::operator[](const int index) const{
 if(index==0) return x;
 if(index==1) return y;
 cout<<"a reserve can only have 2 coordinates, the index must be 0 or 1."<<endl;
 exit(EXIT_FAILURE);
}
float reserve::addReserve(const reserve &r2) const{
 return sqrt(pow(getx()-r2[0],2) + pow(gety()-r2[1],2));
}

float operator+(const reserve &r1, const reserve &r2){
 return r1.addReserve(r2);
}
