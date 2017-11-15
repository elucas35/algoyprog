#include "bintree.h"


//map (leaf's id, leaf of the tree)
map<bintree*,leaf*> mapLeaves ;


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

/*void bintree::readLeaves(){
	vector<reserve> r;
	for(auto it = mapLeaves.cbegin(); it != mapLeaves.cend(); ++it)
  {
	  cout<<"feuille : "<<it->first->root;
    r = (it->second)->getCoordLeaf();
    for(vector<reserve>::iterator k = r.begin(); k!=r.cend(); ++k){
		cout<<"val : x "<<k->getx()<<" y "<<k->gety()<<endl;
	}
  }
  
}*/


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

 float halfpoint;
 halfpoint = (array[left][n]+array[right-1][n]) / 2 ;

 if(left != right - 1 && right != left){ 		
  for(i=left; i<right; i++){
   if(n==0){
    if(array[i][n]<=halfpoint){
     lim = i;
   }
 }else if(n==1) {
  if(array[i][n]>halfpoint){
   lim = i;
 }
}
}
}

t->root=halfpoint; 
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

//searches the closest node to the given coord in the tree
bintree* bintree::searchNode(vector<float> bcoordarray, bintree* t, int depth){

 int dim = 0;
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


//return the list of reserves nearest to the query (which distance is less than the one we had so far)
vector<vector<reserve>> bintree::searchNeighborAreas(vector<int> nearestNeighbors, bintree* indext){
  vector<reserve> rcoordarray;
  vector<vector<reserve>> r;
  unsigned int a =0; int cmpt = 0;
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


/*
void bintree::checkTree(bintree* tree){
  if(tree!=nullptr){
   cout<<" node : "<<tree->root<<endl;
   if( tree->left) bintree::checkTree(tree->left);
   if(tree->right) bintree::checkTree(tree->right);
 }else{
   cout<<"empty tree"<<endl;
 }
 }

 void bintree::checkTreePrev(bintree* indext){
  if(indext->prev){
   cout<< "Le noeud : "<<indext->root<< " a un prev " << indext->prev->root<<endl;
   checkTreePrev(indext->prev);
 }else{
   cout<<"Le noeud : "<<indext->root<<" n'a pas de prev"<<endl;
 }
 }
*/

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
//adds a leaves
void addLeavesToTree(float** reservesArray, bintree* t, int nbreserves){
  int i=0;
  for(i=0; i<nbreserves; i++){
    bintree::add(reservesArray[i], t, 0);
  }
}
