#include <vector>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <string>
#include <iterator>
#include <math.h>
#include <cstring>
#include <map>
#include "queue.h" 
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

static void levels(bintree *);
static bintree *create2DBST(float** array, int depth, int nbreserves, int left = 0, int right = -1, bintree* prevt=nullptr);
static bintree *createhalf2DBST(float** array, int depth, int nbreserves, int left = 0, int right = -1, bintree* prevt=nullptr);
static bintree *createfraction2DBST(float** array, int depth, int nbreserves, int left=0, int right=-1, bintree* prevt=nullptr);
static void add(float* coord, bintree *t, int depth);
static void readBases(ofstream &ofile, bintree* tree, string baseFile, int nbreserves);
static bintree* searchNode(vector<float> bcoordarray, bintree* tree, int depth);
static void openOutputFile(bintree* tree, string outputFile, string baseFile, int nbreserves);
static vector<reserve> searchVectorNode(bintree* indext);
static vector<reserve> readMap();
static vector<float> readMapAll(bintree * indext, vector<float> &c, float dMin);
static bintree * searchNeighbors(bintree* indext);
static bintree * searchFirstLeaf(bintree* indextpp, bool right);
static vector<vector<reserve>> searchNeighborAreas(vector<int> nearestNeighbors, bintree* indext);
//static int searchLeaves(bintree* indextpp, int count);
//static void readBintreeVector(vector<bintree*>);
//static void checkTree(bintree* tree);
//static void checkTreePrev(bintree* indext);
};

// int bintree::searchLeaves(bintree* indextpp, int count){

//   if( indextpp->right== nullptr && indextpp->left == nullptr) {
//     count++;
//   }
//   else{
//     if(indextpp->right) count=searchLeaves(indextpp->right, count);

//     if(indextpp->left) count=searchLeaves(indextpp->left, count);
//   }
//   return count; 
// }


// void bintree::readBintreeVector(vector<bintree*> vec){

//   for (vector<bintree*>::iterator k = vec.begin(); 
//    k != vec.end(); 
//    ++k) 
//   {
//     cout << "feuille : " <<(*k)->root << endl;
//   }

// }

// void bintree::checkTreePrev(bintree* indext){
//  if(indext->prev){
//   cout<< "Le noeud : "<<indext->root<< " a un prev."<<endl;
//   checkTreePrev(indext->prev);
// }else{
//   cout<<"Le noeud : "<<indext->root<<" n'a pas de prev"<<endl;
// }
// }

// void bintree::checkTree(bintree* tree){
//  if(tree!=nullptr){
//   cout<<" node : "<<tree->root<<endl;

//   if( tree->left) bintree::checkTree(tree->left);
//   if(tree->right) bintree::checkTree(tree->right);
// }else{
//   cout<<"empty tree"<<endl;
// }
// }



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

//map (leaf's id, leaf of the tree)
map<bintree*,leaf*> mapLeaves ;

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

//does the partition of a 2D array
int partition2D(float** input, int p, int r, int n)
{
  float pivot = input[r][n];
  int m=0;
  float tmp1=0, tmp2=0;
  if(n==0) {m =1;}
  else if(n==1) {m=0;}

  while ( p < r )
  {
    while ( input[p][n] < pivot )
      p++;

    while ( input[r][n] > pivot )
      r--;


    if ( input[p][n] == input[r][n] )
      p++;

    else if ( p <= r ) {
      tmp1 = input[p][n];
      tmp2 = input[p][m];
      input[p][n] = input[r][n];
      input[p][m] = input[r][m];
      input[r][n] = tmp1;
      input[r][m] = tmp2;
    }
  }
  return r;
}

//recursive implementation of quickSelect algorithm
float quickSelect2D(float** input, int p, int r, int k, int n)
{
  if ( p == r ) return input[p][n]; 
  int j = partition2D(input, p, r, n);
  int length = j - p +1;
  if ( length == k ) return input[j][n];
  else if ( k < length ) return quickSelect2D(input, p, j - 1, k, n);
  else  return quickSelect2D(input, j + 1, r, k - length, n);
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

//adds a leaves
void addLeavesToTree(float** reservesArray, bintree* t, int nbreserves){
  unsigned int i=0, j=0;
  for(i=0; i<nbreserves; i++){
    bintree::add(reservesArray[i], t, 0);
  }
}

//sets the limits of an area
void attributeLimitLeaves(){
  
  vector<reserve> reserves;

  for(auto it = mapLeaves.cbegin(); it != mapLeaves.cend(); ++it)
  {
    reserves = (it->second)->getCoordLeaf();
    (it->second)->setMinx(min2D(reserves, 0));
    (it->second)->setMiny(min2D(reserves, 1));
    (it->second)->setMaxx(max2D(reserves, 0));
    (it->second)->setMaxy(max2D(reserves, 1));
  }
}

//sets the coordinates of a reserve and adds it to the vector of reserves
void setReserveWithCoord(vector<float> coordarray, vector<reserve>& reserves){

  reserve r;
  r.setx(coordarray[0]);
  r.sety(coordarray[1]);
  reserves.push_back(r);

}

//creates an array containing properly all the coordinates of the reserves
float** setReservesArray(vector<reserve>& reserves){

  int index = 0;
  float ** reservesArray = new float *[reserves.size()];
  for(unsigned int i = 0; i <reserves.size(); i++)
   { reservesArray[i] = new float[2];
   }

   for (vector<reserve>::iterator k = reserves.begin(); 
     k != reserves.end(); 
     ++k) 
   {
    index = distance(reserves.begin(), k);
    reservesArray[index][0] = k->getx();
    reservesArray[index][1] = k->gety();
  }
  return reservesArray;
}

//reads an 2D array
void read(float **array, int nbreserves){

  unsigned int i=0, j=0;
  for(i=0; i<nbreserves; i++){
    for(j=0;j<(sizeof(*array)/sizeof(**array)) ;j++){
      cout << array[i][j] << " ";
    }
    cout << endl;
  }
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

//adds coordinates to the correspondant leaf of the map. The id of the leaf is the adress of the correspondent root in the tree
void bintree::add(float* coord, bintree *t, int depth){

 int dim = 0;
 if(depth%2 != 0) { dim = 1; }

 if((t->left)!=nullptr || (t->right)!=nullptr){

  if((coord[dim]<=(t->root))){

   if((t->left)!=nullptr){
    depth++;
    add(coord, t->left, depth);
  }else{
    depth++;
    add(coord, t->right, depth);
  }
}else{

 if((t->right)!=nullptr){
  depth++;
  add(coord, t->right, depth);
}else{
 depth++;
 add(coord, t->left, depth);
}
}
}else if((t->left)==nullptr && (t->right)==nullptr){

 map<bintree*, leaf*>::iterator found = mapLeaves.find(t);

 if(found==mapLeaves.end()){
  leaf* l = new leaf ;
  l->addToLeaf(coord);
  mapLeaves.insert ( pair<bintree*,leaf*>(t, l ));

}else{
  ((found->second))->addToLeaf(coord);
}
}
}

//recursive implementation of quickSort algorithm
void quickSort2D(float **arr, int left, int right, int n) {
  
  int i = left, j = right;
  float tmp1 = 0.0, tmp2 =0.0;
  unsigned int m=0;
  float pivot = arr[(left + right) / 2][n];
  if(n==0) {m =1;}
  else if(n==1) {m=0;}

      /* partition */
  while (i <= j) {
    while (arr[i][n] < pivot)
      i++;
    while (arr[j][n] > pivot)
      j--;
    if (i <= j) {
      tmp1 = arr[i][n];
      tmp2 = arr[i][m];
      arr[i][n] = arr[j][n];
      arr[i][m] = arr[j][m];
      arr[j][n] = tmp1;
      arr[j][m] = tmp2;
      i++;
      j--;
    }
  };

      /* recursion */
  if (left < j)
    quickSort2D(arr, left, j, n);
  if (i < right)
    quickSort2D(arr, i, right, n);
}

//creates a 2D tree using median heuristic spliting
bintree *bintree::create2DBST(float **array, int depth, int nbreserves, int left, int right, bintree* prevt)
{
  int n=0;
  depth++;
  bintree *t = new bintree;
  bintree *tl;
  bintree *tr;

  if(depth%2 != 0) { n = 1; }
  if (right == -1) {right=nbreserves;}

  quickSort2D(array, left, right -1, n);

  if (left == right) { 
   return nullptr; 
 }
 if( left == right - 1 ) {
   return new bintree(array[left][n], prevt);
 }


 int med = (left + right) / 2;
 t->root = array[med][n];
 t->prev = prevt;
 tl = create2DBST(array, depth, med, left, med, t);
 tr = create2DBST(array, depth, right, med + 1, right, t);
 t->left = tl;
 t->right = tr;

 return t;
}

//creates a 2D tree using half heuristic spliting
bintree *bintree::createhalf2DBST(float** array, int depth, int nbreserves, int left, int right, bintree* prevt){ 
  int n=0, i=0, lim=0;
  depth++;
  bintree *t = new bintree;
  bintree *tl;
  bintree *tr;

  if(depth%2 != 0) {n = 1;}

  if (right == -1) { right = nbreserves;}
  quickSort2D(array, left, right-1, n);
  if (left == right) { 
   return nullptr; 
 }
 if (left == right - 1 ) {
   return new bintree(array[left][n], prevt);
 }

 float halfpoint[2];
 halfpoint[0] = (array[left][0]+array[right-1][0]) / 2 ;
 halfpoint[1] = (array[left][1]+array[right-1][1]) / 2 ;

 if(left != right - 1 && right != left){ 		
  for(i=left; i<right; i++){
   if(n==0){
    if(array[i][n]<=halfpoint[n]){
     lim = i;
   }
 }else if(n==1) {
  if(array[i][n]>halfpoint[n]){
   lim = i;
 }
}
}
}

t->root=halfpoint[n]; 
t->prev = prevt;
tl = createhalf2DBST(array, depth, lim, left, lim, t);
tr = createhalf2DBST(array, depth, right, lim +1, right, t);
t->left = tl;
t->right = tr;
return t;
}

//creates a 2D tree using fraction heuristic spliting
bintree *bintree::createfraction2DBST(float** array, int depth, int nbreserves, int left, int right, bintree* prevt){
  int n=0,i=0,lim=0;
  depth++;
  bintree *t = new bintree;
  bintree *tl;
  bintree *tr;

  if(depth%2 != 0) {n = 1;}

  if (right == -1) { right = nbreserves;}
  quickSort2D(array, left, right-1, n);
  if (left == right) { 
   return nullptr; 
 }
 if (left == right - 1 ) {
   return new bintree(array[left][n], prevt);
 }

 float fractionpoint;
 if(left != right - 1 && right != left){ 
   if(right-left>=3){
    for(i=left;i<(left+3);i++){
     fractionpoint+=array[i][n];
     lim=i;
   }
   fractionpoint/=3;
 }else{
   for(i=left;i<right;i++){
    fractionpoint+=array[i][n];
    lim=i;
  }
  fractionpoint/=(right-left);
}

}
t->root=fractionpoint; 
t->prev = prevt;
tl = createfraction2DBST(array, depth, lim, left, lim, t);
tr = createfraction2DBST(array, depth, right, lim +1, right, t);
t->left = tl;
t->right = tr;
return t;
}

//reads the levels of a binary tree
void bintree::levels(bintree *tree)
{
 cout << "levels: ";
 char const *comma = "";
 queue<bintree const *> q;

 for (q.enqueue(tree); !q.empty(); ) {
  bintree const *sub = q.dequeue();

  if (sub->left)
    q.enqueue(sub->left);

  if (sub->right)
    q.enqueue(sub->right);

  std::cout << comma;
  std::cout << sub->root;
  comma = ", ";
}
cout << endl;
}

//searches the closest node to the given coord in the tree
bintree* bintree::searchNode(vector<float> bcoordarray, bintree* t, int depth){

 int dim = 0;
 bool seen = false;
 bintree* e;
 if((depth % 2)!=0){dim=1;}

 if(bcoordarray[dim]<=(t->root)){

   if((t->left)!=nullptr){
                //cout<<"coord : "<< dim<<" "<<bcoordarray[dim]<< " cote gauche"<<endl;
    return searchNode(bcoordarray, (t->left), ++depth);
  }else if((t->right)!=nullptr){
                //cout<<"coord : "<< dim<<" "<<bcoordarray[dim]<< " cote droit"<<endl;
    return searchNode(bcoordarray, (t->right), ++depth);
  }

}else {

 if((t->right)!=nullptr){
                //cout<<"coord : "<< dim<<" "<<bcoordarray[dim]<< " cote droit"<<endl;
  return searchNode(bcoordarray, (t->right), ++depth);
}else if((t->left)!=nullptr){
                //cout<<"coord : "<< dim<<" "<<bcoordarray[dim]<< " cote gauche"<<endl;
  return searchNode(bcoordarray, (t->left), ++depth);
}
}

return t;
}

//searches the first leaf
bintree* bintree::searchFirstLeaf(bintree* indextpp, bool right){

  if(right) {

    if( (indextpp->right) && (indextpp->right->right == nullptr && indextpp->right->left == nullptr)) {

      return indextpp->right;
    }
    else if(indextpp->right) {

     return searchFirstLeaf(indextpp->right, false);
   }else{

    return searchFirstLeaf(indextpp, false);
  }
} else {

  if((indextpp->left) && (indextpp->left->right == nullptr && indextpp->left->left == nullptr)) {

   return indextpp->left;
 }
 else if(indextpp->left){

  return searchFirstLeaf(indextpp->left, true);
}else{

 return searchFirstLeaf(indextpp, true);
}
}  
}

//searches the neighbors areas
bintree* bintree::searchNeighbors(bintree* indext){
 bintree* indextVect;

 if( (indext->prev->left) && ((indext->prev)->left != indext)) 
 {
   if(((indext->prev)->left)->left == nullptr && ((indext->prev)->left)->right == nullptr)
   {
     indextVect=indext->prev->left;

   }else{
     indextVect=bintree::searchFirstLeaf(indext->prev->left, true);
   }
 }
 else if( (indext->prev->right) && (indext->prev)->right != indext)
 {  
   if( (indext->prev->right)->left == nullptr && (indext->prev->right)->right == nullptr)
   { 
     indextVect=indext->prev->right;
   }
   else{
     indextVect=bintree::searchFirstLeaf(indext->prev->right, false);
   }
 }
 else {
   indext = indext->prev;
   bool found=false;
   while(!found){
    if((indext->prev) && ((indext->prev->left != indext) || (indext->prev->right != indext))){

      if( (indext->prev->left) && (indext->prev->left != indext)){

       indextVect=searchFirstLeaf(indext->prev, false);

       found=true;
     }else if( (indext->prev->right) && (indext->prev->right != indext)){

       indextVect= searchFirstLeaf(indext->prev, true);
       found=true;
     }
   }
   indext = indext->prev;
 }
}

return indextVect;
}

//return the list of reserves nearest to the query (which distance is less than the one we had so far)
vector<vector<reserve>> bintree::searchNeighborAreas(vector<int> nearestNeighbors, bintree* indext){
  vector<reserve> rcoordarray;
  vector<vector<reserve>> r;
  unsigned int a =0, cmpt = 0;
  for(auto it = mapLeaves.cbegin(); it != mapLeaves.cend(); ++it)
  {
    if((it->first)!=indext){
     cmpt ++;
     for(a=0;a<nearestNeighbors.size();a++){
      if(cmpt == nearestNeighbors[a]){
        rcoordarray = it->second ->getCoordLeaf() ;
        r.push_back(rcoordarray);
        cout << "Number of reserves lesser than dMin " << r.size() << endl;
      } 
    }

  }


}



return r;
}

//given the adress of the node returned by searchNode, searchVectorNode search into mapLeaves 
//the correspondent leaf of the tree and the vector of reserves into that leaf
vector<reserve> bintree::searchVectorNode(bintree* indext){
 vector<reserve> rcoordarray;
        //cout << "indext "<< indext->root << endl;
 for(auto it = mapLeaves.cbegin(); it != mapLeaves.cend(); ++it)
 {
  if((it->first)==indext){
      					
   rcoordarray =(it->second)->getCoordLeaf();

   break;
 }
}

return rcoordarray;
}

//gets the last leaf in the map
vector<reserve> bintree::readMap(){
  vector<reserve> reserves;
  for(auto it = mapLeaves.cbegin(); it != mapLeaves.cend(); ++it)
  {
    reserves =(it->second)->getCoordLeaf();
  }
  return reserves; 
}

//calculates the distances to all othes areas (areas different from the one where the query belongs)
vector<float> bintree::readMapAll(bintree * indext, vector<float> &c, float dMin){
  vector<float> distToAreas;
  vector<float> vecTemp = {0.0,0.0};
  float temp = 0.0;
  int cmpt = 0;

  for(auto it = mapLeaves.cbegin(); it != mapLeaves.cend(); ++it)
  {
    cmpt++;
    if(it->first != indext){
                //cout << "feuille " << it->first-> root << endl;
                //cout << "in feuille : " << it->first->root << " of index " << cmpt << endl;
     if(c[0]<= (it->second)->getMinx() && c[1]<= (it->second)->getMiny() ){
                  //cout << "vectTempx "<< vecTemp[0] << endl;
                  //cout << "vectTempy "<<vecTemp[1] << endl;
      vecTemp[0] = ( (it->second)->getMinx() );
      vecTemp[1] = ((it->second)->getMiny());
      temp = dist(c, vecTemp);

                  // cout << "min x is "<< (it->second)->getMinx() << endl;
                  // cout << "min y is "<< (it->second)->getMiny() << endl;
                  // cout << "vectTempx "<< vecTemp[0] << endl;
                  // cout << "vectTempy "<<vecTemp[1] << endl;
                  // cout << "consulta x "<< c[0] << endl;
                  // cout << "consulta y "<<c[1] << endl;
                  // cout << "dist is "<< dist(c, vecTemp) << endl;
                  //cout << " cas 2 et 9 distances est : " << temp << endl;
    }else if (c[0]<= (it->second)->getMinx() &&  c[1]>= (it->second)->getMaxy() ){
      vecTemp[0] = ((it->second)->getMinx());
      vecTemp[1] = ((it->second)->getMaxy());
      temp = dist(c, vecTemp);
                  //cout << " cas 3 et 12 distances est : " << temp << endl;
    }else if(c[0]>= (it->second)->getMaxx() &&  c[1]<= (it->second)->getMiny() ){
      vecTemp[0] = ((it->second)->getMaxx());
      vecTemp[1] = ((it->second)->getMiny());
      temp = dist(c, vecTemp);
                  //cout << " cas 5 et 8 distances est : " << temp << endl;
    }else if(c[0]>= (it->second)->getMaxx() &&  c[1]>= (it->second)->getMaxy() ){
      vecTemp[0] = ((it->second)->getMaxx());
      vecTemp[1] = ((it->second)->getMaxy());
      temp = dist(c, vecTemp);
                  //cout << " cas 6 et 11 distances est : " << temp << endl;
    }else if (c[0]<= (it->second)->getMinx() &&  (c[1]<= (it->second)->getMaxy() && c[1]>= (it->second)->getMiny()) ){
      vecTemp[0] = ((it->second)->getMinx());
      vecTemp[1] = (c[1]);
      temp = dist(c, vecTemp);
                  //cout << " cas 1 distances est : " << temp << endl;
    }else if (c[0]>= (it->second)->getMaxx() &&  (c[1]<= (it->second)->getMaxy() && c[1]>= (it->second)->getMiny()) ){
      vecTemp[0] = ((it->second)->getMaxx());
      vecTemp[1] = (c[1]);
      temp = dist(c, vecTemp);
                  //cout << " cas 4 distances est : " << temp << endl;
    }else if (c[1]<= (it->second)->getMiny() &&  (c[0]<= (it->second)->getMaxx() && c[1]>= (it->second)->getMinx()) ){
      vecTemp[0] = (c[0]);
      vecTemp[1] = ((it->second)->getMiny());
      temp = dist(c, vecTemp);
                  //cout << " cas 7 distances est : " << temp << endl;
    }else if (c[1]>= (it->second)->getMaxy() &&  (c[0]<= (it->second)->getMaxx() && c[1]>= (it->second)->getMinx()) ){
      vecTemp[0] = (c[0]);
      vecTemp[1] = ((it->second)->getMaxy());
      temp = dist(c, vecTemp);
                  //cout << " cas 10 distances est : " << temp << endl;
    }else {
      cout << "ATTENTION RENTRE NUL PART point : " << c[0] << " " << c[1] << endl;
    }

    distToAreas.push_back(temp);
  }

}
return distToAreas; 
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

//reads the coordinates provided
vector<float> readCoord(string line, bool &seen){

  char *cstr = new char[line.length() + 1];
  unsigned int i=0, l=0;
  float coord=0.0;
  string s="", t="";
  vector< string > arr;
  vector< float > coordarray;
  bool error = false;


  strcpy(cstr, line.c_str());

  for(i=0;i<line.size();i++){
    if(cstr[i] != ' ') {
      if((isdigit(cstr[i]) || cstr[i] == '.') || cstr[i] == '-') s += cstr[i];
      else  if (seen) {
        cout << "Please replace the character " << cstr[i] << "." <<endl;
      }
      else {
        seen = true;
        cout << "ERROR : The input must contain only characters '.' and '-' and/or numbers." << endl;
        cout << "Please replace the character " << cstr[i] << ". And check the rest of the file."  << endl;
        error = true;
      }
    }
    else {

      t+=s;
      arr.push_back(s);
      s="";
    }

  }
  arr.push_back(s);
  for (l=0; l<arr.size(); l++){

    if( arr[l] != "") {
      coord = stof(arr[l]);
      coordarray.push_back(coord);
    }
  }

  delete [] cstr;
  if(error) {
    cout << "And re-run the program" << endl; exit (EXIT_FAILURE);
  }
  return coordarray;
}

//reads the coordinates provided
vector<float> readCoord(string line, bool &seen, string &filename){

  char *cstr = new char[line.length() + 1];
  unsigned int i=0, l=0;
  float coord = 0.0;
  string s="", t="";
  vector< string > arr;
  vector< float > coordarray;
  bool error = false;


  strcpy(cstr, line.c_str());

  for(i=0;i<line.size();i++){
    if(cstr[i] != ' ') {
      if((isdigit(cstr[i]) || cstr[i] == '.') || cstr[i] == '-') s += cstr[i];
      else  if (seen) {
        cout << "Please replace the character " << cstr[i] << "." <<endl;
      }
      else {
        seen = true;
        cout << "ERROR : The file " << filename << " must contain only the dot character which is '.' and/or numbers." << endl;
        cout << "Please replace the character " << cstr[i] << ". And check the rest of the file."  << endl;
        error = true;
      }
    }
    else {

      t+=s;
      arr.push_back(s);
      s="";
    }

  }
  arr.push_back(s);
  for (l=0; l<arr.size(); l++){

    if( arr[l] != "") {
      coord = stof(arr[l]);
      coordarray.push_back(coord);
    }
  }

  delete [] cstr;
  if(error) {
    cout << "And re-run the program" << endl; exit (EXIT_FAILURE);
  }
  return coordarray;
}

//parses each line through readcoord and returns the number of reserves
int readReserves(vector<reserve>& reserves, string reserveFile){
  string line="", d="";
  int nbreserves = 0;
  unsigned int dimension = 0;
  bool seen = false;
  vector<float> c, output, cd;

  ifstream file;
  file.open(reserveFile);
  if (file.is_open())
  {
    getline (file,d);
    if(d == "" || d== " "){cout << "ERROR : Empty file." << endl; exit (EXIT_FAILURE); }
    else {
      dimension = stoi(d);
      while(getline(file, line)){

        nbreserves++;
        c = readCoord(line, seen, reserveFile);
        if(c.size() != dimension) {
          cout << "ERROR : One reserve has coordinates of wrong dimension, dimension must be "<< dimension <<  " . Please correct it and re-run the program." << endl; exit (EXIT_FAILURE);
        }
        if(!c.size()) {
          cout << "ERROR : the file "<< reserveFile << " must not contain any empty line, please fill the empty line and re-run the program." << endl; exit (EXIT_FAILURE);
        }

        setReserveWithCoord(c, reserves);
      }
      file.close();
    }
  }

  else cout << "Error when opening reserves files " << endl;;
  return nbreserves;
}

//writes the nearest reserve found into an output file
void writeOutputFile(vector<float> output, ofstream &file){
 unsigned int a=0;
 for(a=0; a < output.size(); a++){
  file << output[a] << " ";
}
file << endl;
}

//case 1 : input file + output file 
void bintree::readBases(ofstream &ofile, bintree* tree, string baseFile, int nbreserves){

  string bline="";
  vector<float> output, cd, c, distToAreas;
  vector<int> nearestNeighbors;
  ifstream bfile ;
  bfile.open(baseFile);
  bool seen =false;
  unsigned int dim =2;
  int nblignes=0;
  float dMin  = 0.0;

  if (bfile.is_open())
  {
            //AA
    while(getline(bfile, bline)){
      if(bline == "" || bline== " "){cout << "ERROR : Empty file." << endl; exit (EXIT_FAILURE); }


      c = readCoord(bline, seen, baseFile);
      if(c.size()!=dim){
       cout << "ERROR : The base file has bases with a different dimension than the reserves, dimension must be " << dim << endl; exit (EXIT_FAILURE);
     }
     if(!c.size()) {
      ofile.close();
      ofile.open("output.txt", std::ofstream::out | std::ofstream::trunc);
      ofile << "ERROR : One line is empty, please fill the empty line and re-run the program." << endl;
      cout << "ERROR : the file must not contain any empty line, please fill the empty line and re-run the program." << endl; exit (EXIT_FAILURE);
    }
    if(nbreserves==1||nbreserves==2) {

      vector<reserve> reserves;
      reserves = bintree::readMap();
      cd = evaluateDist(c, reserves);
      dMin = min(cd);
      output = findNearestReserve(min(cd), reserves);
      writeOutputFile(output, ofile);

    }
    else{

      cout << endl << endl << endl;
      bintree* indext = bintree::searchNode(c, tree, 0);
      cout << "feuille initiale " << indext-> root << endl;
      bintree* indextNeighbor= bintree::searchNeighbors(indext);
      vector<reserve> reserves = bintree::searchVectorNode(indext);
      cd = evaluateDist(c, reserves);
      cout << "initial reserves vector size : " << reserves.size() << endl;
      dMin = min(cd);
      cout << "Distances in its own area : " << endl;
      readVecFloat(cd);
      cout << "Min dist so far is : " << cd[dMin] << endl;
      distToAreas = bintree::readMapAll(indext, c, cd[dMin]);
      cout << "Distances to others areas : " << endl;
      readVecFloat(distToAreas);
      nearestNeighbors = lesser(distToAreas, cd[dMin]);
      cout << "reserves in the radius "<< nearestNeighbors.size() << endl;
      vector<vector<reserve>> neighborReserves = bintree::searchNeighborAreas(nearestNeighbors, indext);
      for(vector<vector<reserve>>::iterator k = neighborReserves.begin(); 
       k != neighborReserves.end(); 
       ++k){
        reserves.insert(reserves.end(), k->begin(), k->end());
    }
    
      cd = evaluateDist(c, reserves);
      output = findNearestReserve(min(cd), reserves);
      writeOutputFile(output, ofile);
  }
}
if(!c.size()) {
 cout << "ERROR : Empty file." << endl; exit (EXIT_FAILURE);
}
bfile.close();
}

else cout << "Error when opening bases files " << endl;;
}

//writes the output to the standard output (the screen)
void writeOutput(vector<float> output){

  unsigned int a=0;
  cout <<"The output is : " << endl;
  for(a=0; a < output.size(); a++){
    cout << output[a] << " ";
  }
  cout << endl;
  cout <<"End of output." << endl;
}

//case 2 : input file provided + standard output 
void readBasesNoOutput(bintree* tree, string baseFile, int nbreserves){
  cout << "The nearest reserves to the bases are : " << endl;
  string bline="";
  vector<float> output, cd, c, distToAreas;
  vector<int> nearestNeighbors;
  ifstream bfile ;
  bfile.open(baseFile);
  bool seen =false;
  unsigned int dim = 2;
  int nblignes=0;
  float dMin  = 0.0;

  if (bfile.is_open())
  {
    while(getline(bfile, bline)){
      c = readCoord(bline, seen, baseFile);
      if(c.size()!=dim){
        cout << "ERROR : The bases have a different dimension then the reserves, dimension must be " << dim << endl; exit (EXIT_FAILURE);
      }

      if(nbreserves==1||nbreserves==2) {


        vector<reserve> reserves = bintree::readMap();
        cd = evaluateDist(c, reserves);
        output = findNearestReserve(min(cd), reserves);
        writeOutput(output);

      }
      else{

      bintree* indext = bintree::searchNode(c, tree, 0);
      //cout << "feuille initiale " << indext-> root << endl;
      bintree* indextNeighbor= bintree::searchNeighbors(indext);
      vector<reserve> reserves = bintree::searchVectorNode(indext);
      cd = evaluateDist(c, reserves);
      cout << "initial reserves vector size : " << reserves.size() << endl;
      dMin = min(cd);
      cout << "Distances in its own area : " << endl;
      readVecFloat(cd);
      cout << "Min dist so far is : " << cd[dMin] << endl;
      distToAreas = bintree::readMapAll(indext, c, cd[dMin]);
      cout << "Distances to others areas : " << endl;
      readVecFloat(distToAreas);
      nearestNeighbors = lesser(distToAreas, cd[dMin]);
      cout << "reserves in the radius "<< nearestNeighbors.size() << endl;
      vector<vector<reserve>> neighborReserves = bintree::searchNeighborAreas(nearestNeighbors, indext);
      for(vector<vector<reserve>>::iterator k = neighborReserves.begin(); 
       k != neighborReserves.end(); 
       ++k){
        reserves.insert(reserves.end(), k->begin(), k->end());
      }
    
      cd = evaluateDist(c, reserves);
      output = findNearestReserve(min(cd), reserves);
      writeOutput(output);
      }
    }
    if(!c.size()) {
     cout << "ERROR : Empty file." << endl; exit (EXIT_FAILURE);
   }
   bfile.close();
 }
 else cout << "Error when opening bases files " << endl;; 
}

//case 1 : call readBases() 
void bintree::openOutputFile(bintree* tree, string outputFile, string baseFile, int nbreserves){

  ofstream ofile;
  ofile.open (outputFile.c_str());
  if (!ofile.is_open())
  {
   cout << "Error when opening output file " << endl;;
 }
 readBases(ofile, tree, baseFile, nbreserves);
 ofile.close();
}

//case 3 et 4 : Bases are provided from standard input
vector<string> readBasesFromInput(){
  string line="";
  vector<string> cinbases;
  cout<<"Please enter bases coordinates, one per line, using 2-dimensions (as for reserves) : "<<endl << "Stop with s."<<endl;
  while(line != "s"){
    getline(cin, line);
    if(!line.size() || line == " " || line == "") {
     cout << "ERROR : Empty line." << endl; exit (EXIT_FAILURE);
   }
   else if (line.size()) cinbases.push_back(line);

 }
 cinbases.pop_back();
 return cinbases;
}

//case 3 : standard input + output file
void readBases(bintree* tree, string outputFile, vector<string>& cinbases, ofstream& myfile, int nbreserves){
  string bline="";
  vector<float> output, cd, c, distToAreas;
  vector<int> nearestNeighbors;
  ifstream bfile;
  bfile.open(outputFile);
  bool seen =false;
  unsigned int a=0, l=0, dim = 2;
  int nblignes=0;
  float dMin  = 0.0;

  if(cinbases.size()){
    for(a=0;a<cinbases.size();a++){
      outputFile = "standard input";
      c = readCoord(cinbases[a], seen, outputFile);
      if(c.size()!=dim){
        cout << "ERROR : The base have a different dimension then the reserves, dimension must be " << dim << endl; exit (EXIT_FAILURE);
      }
      if(!c.size()) {
        cout << "ERROR : the file must not contain any empty line, please fill the empty line and re-run the program." << endl; exit (EXIT_FAILURE);
      }
      if(nbreserves==1||nbreserves==2) {


        vector<reserve> reserves =  bintree::readMap();
        cd = evaluateDist(c, reserves);
        output = findNearestReserve(min(cd), reserves);
        writeOutputFile(output, myfile);

      }
      else{
      bintree* indext = bintree::searchNode(c, tree, 0);
      bintree* indextNeighbor= bintree::searchNeighbors(indext);
      vector<reserve> reserves = bintree::searchVectorNode(indext);
      cd = evaluateDist(c, reserves);
      cout << "initial reserves vector size : " << reserves.size() << endl;
      dMin = min(cd);
      cout << "Distances in its own area : " << endl;
      readVecFloat(cd);
      cout << "Min dist so far is : " << cd[dMin] << endl;
      distToAreas = bintree::readMapAll(indext, c, cd[dMin]);
      cout << "Distances to others areas : " << endl;
      readVecFloat(distToAreas);
      nearestNeighbors = lesser(distToAreas, cd[dMin]);
      cout << "reserves in the radius "<< nearestNeighbors.size() << endl;
      vector<vector<reserve>> neighborReserves = bintree::searchNeighborAreas(nearestNeighbors, indext);
      for(vector<vector<reserve>>::iterator k = neighborReserves.begin(); 
       k != neighborReserves.end(); 
       ++k){
        reserves.insert(reserves.end(), k->begin(), k->end());
      }
    
      cd = evaluateDist(c, reserves);
      output = findNearestReserve(min(cd), reserves);

       for(l=0;l<output.size();l++){
        myfile << output[l] << " ";
      }
      myfile << endl;

    }
  }
}else{
  if (bfile.is_open())
  {
    while(getline(bfile, bline)){
     c = readCoord(bline, seen, outputFile);
     if(!c.size()) {
      cout << "ERROR : the file must not contain any empty line, please fill the empty line and re-run the program." << endl; exit (EXIT_FAILURE);
    }
    if(nbreserves==1||nbreserves==2) {


      vector<reserve> reserves = bintree::readMap();
      cd = evaluateDist(c, reserves);
      output = findNearestReserve(min(cd), reserves);
      writeOutputFile(output, myfile);
    }


    else{
      bintree* indext = bintree::searchNode(c, tree, 0);
      bintree* indextNeighbor= bintree::searchNeighbors(indext);
      vector<reserve> reserves = bintree::searchVectorNode(indext);
      cd = evaluateDist(c, reserves);
      cout << "initial reserves vector size : " << reserves.size() << endl;
      dMin = min(cd);
      cout << "Distances in its own area : " << endl;
      readVecFloat(cd);
      cout << "Min dist so far is : " << cd[dMin] << endl;
      distToAreas = bintree::readMapAll(indext, c, cd[dMin]);
      cout << "Distances to others areas : " << endl;
      readVecFloat(distToAreas);
      nearestNeighbors = lesser(distToAreas, cd[dMin]);
      cout << "reserves in the radius "<< nearestNeighbors.size() << endl;
      vector<vector<reserve>> neighborReserves = bintree::searchNeighborAreas(nearestNeighbors, indext);
      for(vector<vector<reserve>>::iterator k = neighborReserves.begin(); 
       k != neighborReserves.end(); 
       ++k){
        reserves.insert(reserves.end(), k->begin(), k->end());
      }
    
      cd = evaluateDist(c, reserves);
      output = findNearestReserve(min(cd), reserves);
   }
 }

 bfile.close();
}
else cout << "Error when opening bases files " << endl;;
}
}

//case 4 : no input file + no output file
void readBasesNoInputOutput(bintree* tree, vector<string>& cinbases, int nbreserves){
 cout << "The nearest reserves to the bases are : " << endl;
 string bline="";
 vector<float> output, cd, c, distToAreas;
 vector<int> nearestNeighbors;
 bool seen =false;
 unsigned int a=0, dim=2;
 int nblignes=0;
 float dMin  = 0.0;

 if(cinbases.size()){
  for(a=0;a<cinbases.size();a++){
    c = readCoord(cinbases[a], seen);
    if(c.size()!=dim){
      cout << "ERROR : The base have a different dimension then the reserves, dimension must be " << dim << endl; exit (EXIT_FAILURE);
    }
    if(nbreserves==1||nbreserves==2) {


      vector<reserve> reserves = bintree::readMap();
      cd = evaluateDist(c, reserves);
      output = findNearestReserve(min(cd), reserves);
      writeOutput(output);

    }
    else{

      bintree* indext = bintree::searchNode(c, tree, 0);
      bintree* indextNeighbor= bintree::searchNeighbors(indext);
      vector<reserve> reserves = bintree::searchVectorNode(indext);
      cd = evaluateDist(c, reserves);
      cout << "initial reserves vector size : " << reserves.size() << endl;
      dMin = min(cd);
      cout << "Distances in its own area : " << endl;
      readVecFloat(cd);
      cout << "Min dist so far is : " << cd[dMin] << endl;
      distToAreas = bintree::readMapAll(indext, c, cd[dMin]);
      cout << "Distances to others areas : " << endl;
      readVecFloat(distToAreas);
      nearestNeighbors = lesser(distToAreas, cd[dMin]);
      cout << "reserves in the radius "<< nearestNeighbors.size() << endl;
      vector<vector<reserve>> neighborReserves = bintree::searchNeighborAreas(nearestNeighbors, indext);
      for(vector<vector<reserve>>::iterator k = neighborReserves.begin(); 
       k != neighborReserves.end(); 
       ++k){
        reserves.insert(reserves.end(), k->begin(), k->end());
      }
    
      cd = evaluateDist(c, reserves);
      output = findNearestReserve(min(cd), reserves);
      writeOutput(output);
      
    }
  }
}
}


int main (int argc, char* argv[]) {
   bool points=false, input=false, output=false, heuristic=false;
   vector<reserve> reserves;
   float **reservesArray=NULL;
   string basef="", outputf="", splitting="";
   vector <string> cinbases;
   vector <float> outputv;
   bintree *tree ;
   int depth = -1, nbreserves = 0;
   vector<bintree*> areas;

   for (unsigned int i = 1; i < argc; ++i) {
    string arg = argv[i];
    if((arg == "-p")||(arg=="--points")) {
     points=true;
     if (i + 1 < argc) {
      nbreserves = readReserves(reserves, argv[++i]);

      reservesArray = setReservesArray(reserves);
    } else {
     cerr << "-p option requires one argument." << endl;
     exit (EXIT_FAILURE);
   }
  } else if((arg == "-i")||(arg=="--input")) {
   input=true;
   if (i + 1 < argc) { 
    basef=argv[++i];
  }else {
   cerr << "-i option requires one argument." << endl;
   exit (EXIT_FAILURE);
  }
  } else if((arg=="-o")||(arg=="--output")){ 
   output=true;
   if (i + 1 < argc) { 
    outputf=argv[++i];
  }else {
   cerr << "-o option requires one argument." << endl;
   exit (EXIT_FAILURE);
  }
  }else if((arg=="-s")||(arg=="--split")){
   heuristic=true;
   if (i + 1 < argc) { 
    splitting=argv[++i];
  }else {
   cerr << "-s option requires one argument." << endl;
   exit (EXIT_FAILURE);
  }

  }
  }

  if(heuristic){
    if(splitting == "mediana"){
     tree = bintree::create2DBST(reservesArray, depth, nbreserves);
     addLeavesToTree(reservesArray, tree, nbreserves);
        			//bintree::levels(tree); 

   }else if(splitting == "mitad"){
     tree = bintree::createhalf2DBST(reservesArray, depth, nbreserves);
     addLeavesToTree(reservesArray, tree, nbreserves);
        			//bintree::levels(tree); 

   }else if(splitting == "promedio"){
     tree = bintree::createfraction2DBST(reservesArray, depth, nbreserves);
     addLeavesToTree(reservesArray, tree, nbreserves);
        			//bintree::levels(tree); 
   }else{
     cerr << "-s option accepts only mediana, mitad or promedio arguments." << endl;
     exit (EXIT_FAILURE);
   }
  }else{

    tree = bintree::create2DBST(reservesArray, depth, nbreserves);
    addLeavesToTree(reservesArray, tree, nbreserves);
    attributeLimitLeaves();
        		//bintree::levels(tree); 

  }
  if((input) && (output))
  {
   bintree::openOutputFile(tree, outputf, basef, nbreserves);

  }
  else if((input) && (!output))
  { 
    readBasesNoOutput(tree, basef, nbreserves);
  }
  if((!input) && (output))
  {
   cinbases = readBasesFromInput();
   ofstream myfile;
   myfile.open (outputf);
   readBases(tree, outputf, cinbases, myfile, nbreserves);
   myfile.close();
  }
  else if((!input) && (!output))
  { 
   cinbases = readBasesFromInput();
   readBasesNoInputOutput(tree, cinbases, nbreserves);
  }
  if(!points) {cerr << "ERROR : -p or --points is a compulsory option." << endl; exit (EXIT_FAILURE);}
  if (tree != nullptr) {delete tree;} 
  for (unsigned int i = 0;i<reserves.size();i++) { delete[] reservesArray[i]; }
    delete[] reservesArray;

return 0;
}
