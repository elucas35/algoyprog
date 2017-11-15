#ifndef BINTREE_H_INCLUDED
#define BINTREE_H_INCLUDED

#include <map>
#include "leaf.h"
#include "distance.h"
class bintree 
{
private:
  float root;
  bintree *left;
  bintree *right;
  bintree *prev;

public:
 bintree() : left(nullptr), right(nullptr), prev(nullptr){
 }

 bintree(const float &t) : root(t), left(nullptr), right(nullptr), prev(nullptr) {
 }

 bintree(const float &t, bintree* prevt) : root(t), left(nullptr), right(nullptr), prev(prevt) {
 }

 ~bintree() {
  if (left != nullptr)
    delete left;

  if (right != nullptr)
    delete right;

 }

static bintree *create2DBST(float** array, int depth, int nbreserves, int left = 0, int right = -1, bintree* prevt=nullptr);
static bintree *createhalf2DBST(float** array, int depth, int nbreserves, int left = 0, int right = -1, bintree* prevt=nullptr);
static bintree *createfraction2DBST(float** array, int depth, int nbreserves, int left=0, int right=-1, bintree* prevt=nullptr);
static void add(float* coord, bintree *t, int depth);
static bintree* searchNode(vector<float> bcoordarray, bintree* tree, int depth);
static vector<reserve> searchVectorNode(bintree* indext);
static vector<reserve> readMap();
static vector<float> readMapAll(bintree * indext, vector<float> &c, float dMin);
static vector<vector<reserve>> searchNeighborAreas(vector<int> nearestNeighbors, bintree* indext);
//static void checkTree(bintree* tree);
//static void checkTreePrev(bintree* indext);
//static void readLeaves();
};
void quickSort2D(float **arr, int left, int right, int n);
void attributeLimitLeaves();
void addLeavesToTree(float** reservesArray, bintree* t, int nbreserves);
#endif
