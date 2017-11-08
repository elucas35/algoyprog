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
    
    public:
    void setx(float);
    void sety(float);
    float getx();
    float gety();
    int getDim();
  };

  void reserve::setx(float coordx){
    x = coordx;
  }
  void reserve::sety(float coordy){
    y = coordy;
  }
  float reserve::getx(){
    return x;
  }
  float reserve::gety(){
    return y;
  }
  int reserve::getDim(){
    cout << "2" << endl;
    return 2;
}

class bintree {
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

	  static void niveles(bintree *);
	  static bintree *create2DBST(float** array, int depth, int nbreserves, int left = 0, int right = -1, bintree* prevt=nullptr);
	  static bintree *createhalf2DBST(float** array, int depth, int nbreserves, int left = 0, int right = -1);
	  static bintree *createfraction2DBST(float** array, int depth, int nbreserves, int left=0, int right=-1);
	  static void add(float* coord, bintree *t, int depth);
	  static void readBases(ofstream &ofile, bintree* tree, string baseFile);
	  static vector<float> evaluateDist(vector<float> bcoordarray, bintree* tree, vector<reserve> reserves);
	  static bintree* searchNode(vector<float> bcoordarray, bintree* tree, int depth);
	  static void openOutputFile(bintree* tree, string outputFile, string baseFile);
	  static void checkTree(bintree* tree);
	  static void checkTreePrev(bintree* indext);
	  static vector<reserve> searchVectorNode(bintree* indext);
      static vector<bintree *> searchNeighbors(bintree* indext, int nArea);
      static bintree* searchFirstLeaf(bintree* indextpp, bool right);
      static int searchLeaves(bintree* indextpp, int count);
      static void searchOtherNeighbors(bintree* indext, int nArea, vector<bintree*> &indexVector);

    
};

 //leaves of the tree
class leaf{
	private:
		vector< reserve > coordleaf;
	public:
		void addToLeaf(float* coord);
		vector<reserve> getCoordLeaf();
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



void addLeavesToTree(float** reservesArray, bintree* t, int nbreserves){
	 unsigned int i, j;
    for(i=0; i<nbreserves; i++){
		cout<<reservesArray[i]<<" ";
        bintree::add(reservesArray[i], t, 0);
    }
}
 
//map (leaf's id, leaf of the tree)
map<bintree*,leaf*> mapLeaves ;


//calculate the distance between coord of the base and coord of the reserve
float dist (vector<float>& d1, vector<float>& d2){
    unsigned int a;
    float d,sum=0;
    for(a=0; a < d1.size(); a++){
      d=d2[a] - d1[a];
      sum += d*d;
    }
    return sqrt(sum);
}

//return the indice of the smallest distance
int min(vector<float> &d){
    unsigned int a;
    int indice = 0;
    for(a=0; a < d.size(); a++){
     if (d[a] <  d[indice]) {
      indice = a;
     }
    }
    return indice;
}

void setReserveWithCoord(vector<float> coordarray, vector<reserve>& reserves){

    reserve r;
    r.setx(coordarray[0]);
    r.sety(coordarray[1]);
    reserves.push_back(r);

}

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

void read(float **array, int nbreserves){

    unsigned int i, j;
    for(i=0; i<nbreserves; i++){
      for(j=0;j<(sizeof(*array)/sizeof(**array)) ;j++){
        cout << array[i][j] << " ";
      }
      cout << endl;
    }
    cout << "row : " << nbreserves << endl;
    cout << "col : " << sizeof *array / sizeof **array  << endl;
}

vector<float> findNearestReserve(int index, vector<reserve>& reserves){
   // cout<<"index "<<index<<endl;
    vector<float> output={0,0};
    
    vector<reserve>::iterator k = reserves.begin();
    advance (k,index); 
    output[0] = k->getx();
    output[1] = k->gety();
    return output;
}

 void bintree::checkTree(bintree* tree){
	  if(tree!=nullptr){
		  cout<<" noeud : "<<tree->root<<endl;
		 
		  if( tree->left) bintree::checkTree(tree->left);
		  if(tree->right) bintree::checkTree(tree->right);
	  }else{
		  cout<<"arbre vide"<<endl;
	  }
  }
  
  
//Add coordinates to the correspondant leaf of the map. The id of the leaf is the adress of the correspondent root in the tree
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
		  cout<<"Nouvelle feuille ajoutée : "<< t->root<< " "<<t<<" coord : "<<coord[0]<<" , "<<coord[1]<<endl;
		  leaf* l = new leaf ;
		  l->addToLeaf(coord);
		  mapLeaves.insert ( pair<bintree*,leaf*>(t, l ));
		   
	  }else{
		  cout<<"Coord ajoutées à feuille déjà existante : "<< t->root<< " "<<t<<" coord : "<<coord[0]<<" , "<<coord[1]<<endl;
		  ((found->second))->addToLeaf(coord);
	  }
  }
}
//sorting algorithm to place the values of the array in order
void quickSort2D(float **arr, int left, int right, int n) {
      int i = left, j = right;
      float tmp1, tmp2;
      int m;
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

//finds recursively the median values of the binary search tree to build
bintree *
bintree::create2DBST(float **array, int depth, int nbreserves, int left, int right, bintree* prevt)
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
  if(prevt != nullptr) cout  << " prevt " << prevt->root << " profondeur : "<< depth<< endl;
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

  
//   if( left == right - 5 ) {

//     return new bintree(array[(( (right - left) / 2) + left)][n]);
//   }
//   else if(right - left > 5) {
//     cout << "tree"<< endl;
//   int med = (left + right) / 2;
//   t->root = array[med][n];
//   tl = create2DBST(array, depth, med, left, med);
//   tr = create2DBST(array, depth, right, med , right);
//   t->left = tl;
//   t->right = tr;

// }

//finds recursively the half values of the binary search tree to build
bintree *
bintree::createhalf2DBST(float** array, int depth, int nbreserves, int left, int right){ 
  int n = 0, i, lim;
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
	  return new bintree(array[left][n]);
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
    tl = createhalf2DBST(array, depth, lim, left, lim);
    tr = createhalf2DBST(array, depth, right, lim +1, right);
    t->left = tl;
    t->right = tr;
    return t;
}

//finds recursively the half value of a fraction of points of the binary search tree to build
bintree *
bintree::createfraction2DBST(float** array, int depth, int nbreserves, int left, int right){
  int n=0,i,lim;
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
	  return new bintree(array[left][n]);
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
	tl = createfraction2DBST(array, depth, lim, left, lim);
    tr = createfraction2DBST(array, depth, right, lim +1, right);
    t->left = tl;
    t->right = tr;
    return t;
}


void
bintree::niveles(bintree *tree)
{
	//readTreeLeaves();
	cout << "niveles: ";
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

//search the closest node to the given coord in the tree
bintree* bintree::searchNode(vector<float> bcoordarray, bintree* t, int depth){
	
	int dim = 0;
	if((depth % 2)!=0){dim=1;}
	
		
		if(bcoordarray[dim]<=(t->root)){
		 
			  if((t->left)!=nullptr){
				  return searchNode(bcoordarray, (t->left), ++depth);
			  }
			  
		}else {
			
			  if((t->right)!=nullptr){
				  return searchNode(bcoordarray, (t->right), ++depth);
			  }
		}
		
		  return t;
}
void bintree::searchOtherNeighbors( bintree *indext, int nArea, vector<bintree*> &indexVector){
	bool right = false;
	if((indext->prev->left) && (indext->prev->left) != indext){
		right=true;
	}
	if(!right && (indext->prev->prev->right)){
		if(indext->prev->root < indext->prev->prev->right->root){
			
			if((indext->prev->prev->right)->left )cout<< "on va voir l'enfant gauche de "<<indext->prev->prev->right->root<< " qui est " << indext->prev->prev->right->left->root<<endl;
			/*indexVector.push_back(searchFirstLeaf(indext->prev->prev->right, !right);*/
		}else{
			
			if((indext->prev->prev->right)->left ) cout<< "on va voir l'enfant gauche de "<<indext->prev->prev->right->root<< " qui est " << indext->prev->prev->right->left->root;
			if((indext->prev->prev->right)->right ) cout<<" et son enfant droit qui est "<<indext->prev->prev->right->right->root<<endl;
			/*indexVector.push_back(searchFirstLeaf(indext->prev->prev->right, !right);*/
			/*indexVector.push_back(searchFirstLeaf(indext->prev->prev->right, right);*/
		}
	}
	else if(right && (indext->prev->prev->left)){
		if(indext->prev->root < indext->prev->prev->left->root){
			if(indext->prev->prev->left->left) cout<< "on va voir l'enfant gauche de "<<indext->prev->prev->left->root<< " qui est " << indext->prev->prev->left->left->root<<endl;
			/*indexVector.push_back(searchFirstLeaf(indext->prev->prev->left, right);*/
		}else{
			if(indext->prev->prev->left->left) cout<< "on va voir l'enfant gauche de "<<indext->prev->prev->left->root<< " qui est " << indext->prev->prev->left->left->root;
			if(indext->prev->prev->left->right) cout<<" et son enfant droit qui est "<<indext->prev->prev->right->left->root<<endl;
			/*indexVector.push_back(searchFirstLeaf(indext->prev->prev->left, !right);*/
			/*indexVector.push_back(searchFirstLeaf(indext->prev->prev->left, right);*/
		}
	}
}
bintree* bintree::searchFirstLeaf(bintree* indextpp, bool right){
  cout <<"indextpp " <<  indextpp->root <<  right << endl;
  if(right) {
	  cout<<"je vais explorer le coté droit"<<endl;
    if( (indextpp->right) && (indextpp->right->right == nullptr && indextpp->right->left == nullptr)) {
		cout << indextpp->right->root << " est une feuille" << endl; 
		return indextpp->right;
	}
    else if(indextpp->right) {
		cout << indextpp->right->root << " est un noeud" << endl; 
       searchFirstLeaf(indextpp->right, false);
    }else{
		cout<< "Il n'y a pas de coté droit"<<endl;
		searchFirstLeaf(indextpp, false);
	}
   } else {
	  cout<<"je vais explorer le coté gauche"<<endl;
		if((indextpp->left) && (indextpp->left->right == nullptr && indextpp->left->left == nullptr)) {
			cout << indextpp->left->root << " est une feuille" << endl; 
			return indextpp->left;
		}
		else if(indextpp->left){
			cout << indextpp->left->root << " est un noeud" << endl; 
		    searchFirstLeaf(indextpp->left, true);
		}else{
		cout<< "Il n'y a pas de coté gauche"<<endl;
			searchFirstLeaf(indextpp, true);
		}
    }  
}

int bintree::searchLeaves(bintree* indextpp, int count){
	
    if( indextpp->right== nullptr && indextpp->left == nullptr) {
		count++;
	}
    else{
		if(indextpp->right) count=searchLeaves(indextpp->right, count);
       
		if(indextpp->left) count=searchLeaves(indextpp->left, count);
    }
    return count; 
}


//search the neighbors areas
vector<bintree*> bintree::searchNeighbors(bintree* indext, int nArea){
	
  vector<bintree*> indextVect;
  //verifie qu'il y a un frere gauche et qu'il n'est pas indext 
  if( (indext->prev->left) && ((indext->prev)->left != indext)) 
  {
	  //ce frere gauche est une feuille
	  if(((indext->prev)->left)->left == nullptr && ((indext->prev)->left)->right == nullptr)
	  {
		cout << "existe feuille soeur gauche " << endl; 
		indextVect.push_back((indext->prev)->left);
		bintree::searchOtherNeighbors(indext, nArea, indextVect);
		
		//ce frere gauche a des enfants
	  }else{
		indextVect.push_back(bintree::searchFirstLeaf(indext->prev->left, true));
		bintree::searchOtherNeighbors(indext, nArea, indextVect);
	  }
  }
	//verifie qu'il y a un frere droit et qu'il n'est pas indext 
  else if( (indext->prev->right) && (indext->prev)->right != indext)
  {  
	  //ce frere droit est une feuille
	  if( (indext->prev->right)->left == nullptr && (indext->prev->right)->right == nullptr)
		{ 
			cout << "existe feuille soeur droite " << endl;
			indextVect.push_back((indext->prev)->right);
			bintree::searchOtherNeighbors(indext, nArea, indextVect);
		}
		//ce frere droit a des enfants
		else{
		indextVect.push_back(bintree::searchFirstLeaf(indext->prev->right, false));
		bintree::searchOtherNeighbors(indext, nArea, indextVect);
		}
  }
  //n'a pas de frere donc on va explorer les enfants du grand-père du coté où nous n'étions pas
  else {
	  /*cout << "n'a pas de frere" << endl;*/
	  indext = indext->prev;
	  cout<<"racine de mon pere : "<<indext->root<<endl;
	  bool found=false;
	  while(!found){
		  if((indext->prev) && ((indext->prev->left != indext) || (indext->prev->right != indext))){
			  cout<<"il y a un grand pere qui a des chemins différents du mien"<<endl;
				if( (indext->prev->left) && (indext->prev->left != indext)){
					cout<<"racine de mon grandpere : "<<indext->prev->root<<endl;
					indextVect.push_back(searchFirstLeaf(indext->prev, false));
					bintree::searchOtherNeighbors(indext->right, nArea, indextVect);
					found=true;
				}else if( (indext->prev->right) && (indext->prev->right != indext)){
					cout<<"racine de mon grandpere : "<<indext->prev->root<<endl;
					indextVect.push_back(searchFirstLeaf(indext->prev, true));
					bintree::searchOtherNeighbors(indext->left, nArea, indextVect);
					found=true;
				}
			}
			indext = indext->prev;
		}
  }
	
  return indextVect;
}

//given the adress of the node returned by searchNode, searchVectorNode search into mapLeaves 
//the correspondent leaf of the tree and the vector of reserves into that leaf
vector<reserve> bintree::searchVectorNode(bintree* indext){
	vector<reserve> rcoordarray;
	for(auto it = mapLeaves.cbegin(); it != mapLeaves.cend(); ++it)
			{
				if((it->first)==indext){
					
					rcoordarray =(it->second)->getCoordLeaf();
				
					break;
				}
			}
			
			
	return rcoordarray;
}

//return a vector of distances between a base and the reserves
vector<float> bintree::evaluateDist(vector<float> bcoordarray, bintree* tree, vector<reserve> reserves){

    vector<float> rcoordarray = {0,0};
    vector<float> distances;
     for (vector<reserve>::iterator k = reserves.begin(); 
                             k != reserves.end(); 
                             ++k) 
                {
					rcoordarray[0]=k->getx();
					rcoordarray[1]=k->gety();
					//cout<< rcoordarray[0]<< " "<<rcoordarray[1]<<endl;
					
				}
				
	distances.push_back(dist(rcoordarray, bcoordarray));
	//cout<<"distance : "<<dist(rcoordarray, bcoordarray)<<endl;
    return distances;
}
	
vector<float> readCoord(string line, bool &seen){
	
    char *cstr = new char[line.length() + 1];
    unsigned int i, l;
    float coord;
    string s, t;
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
	
	
vector<float> readCoord(string line, bool &seen, string &filename){

    char *cstr = new char[line.length() + 1];
    unsigned int i, l;
    float coord;
    string s, t;
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
    
int readReserves(vector<reserve>& reserves, string reserveFile){
    string line, d;
    int nbreserves = 0;
    unsigned int dimension = 0;
    bool seen =false;
    vector<float> c, output, cd;


    ifstream file;
    file.open(reserveFile);
    if (file.is_open())
    {
      getline (file,d);
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

    else cout << "Error when opening reserves files " << endl;;
    return nbreserves;
}

void writeOutputFile(vector<float> output, ofstream &file){
	unsigned int a;
    for(a=0; a < output.size(); a++){
      file << output[a] << " ";
    }
    file << endl;
}


 void bintree::checkTreePrev(bintree* indext){
	if(indext->prev){
		cout<< "Le noeud : "<<indext->root<< " a un prev."<<endl;
		checkTreePrev(indext->prev);
	}else{
		cout<<"Le noeud : "<<indext->root<<" n'a pas de prev"<<endl;
	}
  }


//When both input and output file are provided
void bintree::readBases(ofstream &ofile, bintree* tree, string baseFile){
	
    string bline;
    vector<float> output, cd, c;
    ifstream bfile ;
    bfile.open(baseFile);
    bool seen =false;
    unsigned int dim =2, nArea =0;
    
	nArea = bintree::searchLeaves(tree,0);
    if (bfile.is_open())
    {
      while(getline(bfile, bline)){
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
		  //recherche du noeud pointant sur la feuille la plus proche de la base c
		  bintree* indext = bintree::searchNode(c, tree, 0);
		  //recuperation du vecteur de coordonnées de réserves stockées sur la feuille
		  vector<reserve> reserves = bintree::searchVectorNode(indext);
		  cd = bintree::evaluateDist(c, tree, reserves);
		  output = findNearestReserve(min(cd), reserves);
		  bintree::checkTreePrev(indext);
		  writeOutputFile(output, ofile);
		  vector<bintree*> indextVector= bintree::searchNeighbors(indext, nArea);
		   for (vector<bintree*>::iterator k = indextVector.begin(); 
                             k != indextVector.end(); 
                             ++k) 
                {
					cout<< (*k)->root <<" "<<endl;
				}
      }

      bfile.close();
    }
    else cout << "Error when opening bases files " << endl;; 
  }
  
  void writeOutput(vector<float> output){
    unsigned int a;
    for(a=0; a < output.size(); a++){
      cout << output[a] << " ";
    }
    cout << endl;
  }
  
//when there is no output file provided
void readBasesNoOutput(bintree* tree, string baseFile){
    cout << "The nearest reserves to the bases are : " << endl;
    string bline;
    vector<float> output, cd, c;
    ifstream bfile ;
    bfile.open(baseFile);
    bool seen =false;
    unsigned int dim = 2;

    if (bfile.is_open())
    {
      while(getline(bfile, bline)){
      c = readCoord(bline, seen, baseFile);
      if(c.size()!=dim){
        cout << "ERROR : The bases have a different dimension then the reserves, dimension must be " << dim << endl; exit (EXIT_FAILURE);
      }
		//recherche du noeud pointant sur la feuille la plus proche de la base c
		bintree* indext = bintree::searchNode(c, tree, 0);
		//recuperation du vecteur de coordonnées de réserves stockées sur la feuille
		vector<reserve> reserves = bintree::searchVectorNode(indext);
		cd = bintree::evaluateDist(c, tree, reserves);
		output = findNearestReserve(min(cd), reserves);
		writeOutput(output);
      }

      bfile.close();
    }
    else cout << "Error when opening bases files " << endl;; 
  }
  
void bintree::openOutputFile(bintree* tree, string outputFile, string baseFile){
    
    ofstream ofile;
    ofile.open (outputFile.c_str());
    if (!ofile.is_open())
    {
  	  cout << "Error when opening output file " << endl;;
    }
    readBases(ofile, tree, baseFile);
    ofile.close();
}

vector<string> readBasesFromInput(){
      string line;
      vector<string> cinbases;
      cout<<"Please enter bases coordinates, one per line, using 2-dimensions (as for reserves) : "<<endl << "Stop with s."<<endl;
      while(line != "s"){
		  getline(cin, line);
		  if (line.size()) cinbases.push_back(line);
		  
      }
      cinbases.pop_back();
    return cinbases;
}

//when there is no input file provided
void readBases(bintree* tree, string outputFile, vector<string>& cinbases, ofstream& myfile){
    string bline;
    vector<float> output, cd, c;
    ifstream bfile ;
    bfile.open(outputFile);
    bool seen =false;
    unsigned int a, l, dim = 2;

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
	  //recherche du noeud pointant sur la feuille la plus proche de la base c
	  bintree* indext = bintree::searchNode(c, tree, 0);
	  //recuperation du vecteur de coordonnées de réserves stockées sur la feuille
	  vector<reserve> reserves = bintree::searchVectorNode(indext);
	  cd = bintree::evaluateDist(c, tree, reserves);
	  output = findNearestReserve(min(cd), reserves);
      
      for(l=0;l<output.size();l++){
        myfile << output[l] << " ";
      }
      myfile << endl;

      }

    }else{
		if (bfile.is_open())
		{
		  while(getline(bfile, bline)){
			c = readCoord(bline, seen, outputFile);
			 if(!c.size()) {
				cout << "ERROR : the file must not contain any empty line, please fill the empty line and re-run the program." << endl; exit (EXIT_FAILURE);
			 }
			  //recherche du noeud pointant sur la feuille la plus proche de la base c
			  bintree* indext = bintree::searchNode(c, tree, 0);
			  //recuperation du vecteur de coordonnées de réserves stockées sur la feuille
			  vector<reserve> reserves = bintree::searchVectorNode(indext);
			  cd = bintree::evaluateDist(c, tree, reserves);
			  output = findNearestReserve(min(cd), reserves);
		  }

		  bfile.close();
		}
    else cout << "Error when opening bases files " << endl;;
  }
}

void readBasesNoInputOutput(bintree* tree, vector<string>& cinbases){
	cout << "The nearest reserves to the bases are : " << endl;
    string bline;
    vector<float> output, cd, c;
    bool seen =false;
    unsigned int a, dim=2;

    if(cinbases.size()){
		for(a=0;a<cinbases.size();a++){
		  c = readCoord(cinbases[a], seen);
      if(c.size()!=dim){
        cout << "ERROR : The base have a different dimension then the reserves, dimension must be " << dim << endl; exit (EXIT_FAILURE);
      }
		  //recherche du noeud pointant sur la feuille la plus proche de la base c
		  bintree* indext = bintree::searchNode(c, tree, 0);
		  //recuperation du vecteur de coordonnées de réserves stockées sur la feuille
		  vector<reserve> reserves = bintree::searchVectorNode(indext);
		  cd = bintree::evaluateDist(c, tree, reserves);
		  output = findNearestReserve(min(cd), reserves);
		  writeOutput(output);
		}
	}
  }


int main (int argc, char* argv[]) {
  	bool points=false, input=false, output=false, heuristic=false;
  	vector<reserve> reserves;
    float **reservesArray=NULL;
  	string basef, outputf, splitting;
    vector <string> cinbases;
    vector <float> outputv;
    bintree *tree ;
    int depth = -1;
  	int nbreserves;
  	
  	for (int i = 1; i < argc; ++i) {
          string arg = argv[i];
  		if((arg == "-p")||(arg=="--points")) {
  			points=true;
  			if (i + 1 < argc) {
  				nbreserves = readReserves(reserves, argv[++i]);
				reservesArray = setReservesArray(reserves);
				//read(reservesArray, nbreserves);
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
			bintree::niveles(tree); //Recorrido: por niveles
			
		}else if(splitting == "mitad"){
			bintree *tree = bintree::createhalf2DBST(reservesArray, depth, nbreserves);
      addLeavesToTree(reservesArray, tree, nbreserves);
			bintree::niveles(tree); //Recorrido: por niveles
			
		}else if(splitting == "promedio"){
			bintree *tree = bintree::createfraction2DBST(reservesArray, depth, nbreserves);
      addLeavesToTree(reservesArray, tree, nbreserves);
			bintree::niveles(tree); //Recorrido: por niveles
			
		}else{
			cerr << "-s option accepts only mediana, mitad or promedio arguments." << endl;
            exit (EXIT_FAILURE);
		}
	}else{
		
		tree = bintree::create2DBST(reservesArray, depth, nbreserves);
		addLeavesToTree(reservesArray, tree, nbreserves);
		bintree::niveles(tree); //Recorrido: por niveles
		
	}
	if((input) && (output))
  	 {
  	  bintree::openOutputFile(tree, outputf, basef);
  	 }
  	 else if((input) && (!output))
  	 { 
      readBasesNoOutput(tree, basef);
  	 }
  	 if((!input) && (output))
  	 {
       cinbases = readBasesFromInput();
       ofstream myfile;
       myfile.open (outputf);
       readBases(tree, outputf, cinbases, myfile);
       myfile.close();
       cout<<"The answers are in the output file you provided."<<endl;
    }
    else if((!input) && (!output))
  	 { 
       cinbases = readBasesFromInput();
       readBasesNoInputOutput(tree, cinbases);
	 }
    if(!points) {cerr << "ERROR : -p or --points is a compulsory option." << endl; exit (EXIT_FAILURE);}
  	 if (tree != nullptr) {delete tree;} //Destrucción
    for (unsigned int i = 0;i<reserves.size();i++) { delete[] reservesArray[i]; }
    delete[] reservesArray;
    
    return 0;
}
