#ifndef DISTANCE_H_INCLUDED
#define DISTANCE_H_INCLUDED
#include <vector>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <string>
#include <iterator>
#include <math.h>
#include <cstring>
#include "reserve.h"
using namespace std;

float min2D(vector<reserve> reserves, int dim);
float max2D(vector<reserve> reserves, int dim);
void readVecFloat(vector<float> distToAreas);
float dist (vector<float>& d1, vector<float>& d2);
vector<int> lesser(vector<float> &d, float dMin);
int min(vector<float> &d);
vector<float> findNearestReserve(int index, vector<reserve>& reserves);
vector<float> evaluateDist(vector<float> bcoordarray, vector<reserve> reserves);
#endif
