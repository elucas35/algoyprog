#include "distance.h"
//returns the minimum of an 2D array
float min2D(vector<reserve> reserves, int dim){
  float indice = 0.0;
  float coord = 0.0;
  vector<reserve>::iterator k = reserves.begin();

  if (dim == 0 ) indice = k->getx();
  else indice = k->gety();

  for (vector<reserve>::iterator k = reserves.begin(); 
   k != reserves.end(); 
   ++k) 
  {
    if (dim == 0 ) coord = k->getx();
    else coord = k->gety();
    if (coord <  indice) {
      indice = coord;
    }
  }
  return indice;
}


//returns the maximum of an 2D array
float max2D(vector<reserve> reserves, int dim){
  float indice = 0.0;
  float coord = 0.0;
  vector<reserve>::iterator k = reserves.begin();

  if (dim == 0 ) indice = k->getx();
  else indice = k->gety();

  for (vector<reserve>::iterator k = reserves.begin(); 
   k != reserves.end(); 
   ++k) 
  {
    if (dim == 0 ) coord = k->getx();
    else coord = k->gety();
    if (coord >  indice) {
      indice = coord;
    }
  }
  return indice;
}
//reads a vector of float
void readVecFloat(vector<float> distToAreas){
  unsigned int a=0;
  for(a=0;a<distToAreas.size();a++){
    cout << "Distances are : " << distToAreas[a] << endl;
  }
}
//calculate the distance between 2 pts : a query coordinates and a reserve
float dist (vector<float>& d1, vector<float>& d2){
  unsigned int a=0;
  float d=0.0,sum=0.0;
  for(a=0; a < d1.size(); a++){
    d=d2[a] - d1[a];
    sum += d*d;
  }
  return sqrt(sum);
}
//return the indexes of the values smaller than a specific distance
vector<int> lesser(vector<float> &d, float dMin){
  unsigned int a=0;
  vector<int> indexes;
  for(a=0; a < d.size(); a++){
   if (d[a] <=  dMin) {
    indexes.push_back(d[a]);
  }
}
return indexes;
}

//finds the reserve of index given in a reserve vector, index is of smallest distance found
vector<float> findNearestReserve(int index, vector<reserve>& reserves){

  vector<float> output={0.0,0.0};

  vector<reserve>::iterator k = reserves.begin();
  advance (k,index); 
  output[0] = k->getx();
  output[1] = k->gety();
  return output;
}
//return the index of the smallest distance
int min(vector<float> &d){
  unsigned int a=0;
  int indice = 0;
  for(a=0; a < d.size(); a++){
   if (d[a] <  d[indice]) {
    indice = a;
  }
}
return indice;
}


//returns a vector of distances between a base and the reserves
vector<float> evaluateDist(vector<float> bcoordarray, vector<reserve> reserves){

  vector<float> rcoordarray = {0,0};
  vector<float> distances;
  for (vector<reserve>::iterator k = reserves.begin(); 
   k != reserves.end(); 
   ++k) 
  {
    rcoordarray[0]=k->getx();
    rcoordarray[1]=k->gety();
    distances.push_back(dist(rcoordarray, bcoordarray));
  }

  return distances;
}
