#include <vector>
#include <cstdlib>
#include <iostream>
#include "queue.h" 
#include <fstream>
#include <string>
#include <iterator>
#include <math.h>
#include <cstring>
#include <map>
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

  class bintree {
	private:
	  float root;
	  bintree *left;
	  bintree *right;

	public:
	  bintree() : left(nullptr), right(nullptr){
	  }

	  bintree(const float &t) : root(t), left(nullptr), right(nullptr) {
	  }

	  ~bintree() {
		if (left != nullptr)
		  delete left;

		if (right != nullptr)
		  delete right;
	  }

	  static void niveles(bintree *);
	  static bintree *create2DBST(float** array, int depth, int nbreserves, int left = 0, int right = -1);
	  static bintree *createhalf2DBST(float** array, int depth, int nbreserves, int left = 0, int right = -1);
	  static bintree *createfraction2DBST(float** array, int depth, int nbreserves, int left=0, int right=-1);
	  bool add(float* coord, float* index);
	  static void readBases(ofstream &ofile, bintree* tree, string baseFile);
	  static vector<float> evaluateDist(vector<float> bcoordarray, bintree* tree, vector<reserve> reserves);
	  static float* searchNode(vector<float> bcoordarray, bintree* tree, int depth);
	  static void openOutputFile(bintree* tree, string outputFile, string baseFile);
};

	map<float*,leaf*> mapLeaves ;
	
	
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
  
  void read(float **array)
  {

    unsigned int i, j;
    for(i=0; i<(sizeof(array)/sizeof(array[0])); i++){
      for(j=0;j<(sizeof(*array)/sizeof(**array)) ;j++){
        cout << array[i][j] << " ";
      }
      cout << endl;
    }
    int r = sizeof array / sizeof *array;
    cout << "row : " << r << endl;
    cout << "col : " << sizeof *array / sizeof **array  << endl;

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
  
  
void readTreeLeaves(){
	vector<reserve> v;
	   for(auto it = mapLeaves.cbegin(); it != mapLeaves.cend(); ++it)
			{
				cout << "leaf : "<< *(it->first) << " vector : " ;
				 v =((it->second)->getCoordLeaf());
			
				for(auto it2=v.cbegin();it2!=v.cend();++it2){
					reserve r = *it2;
					float x = r.getx();
					float y = r.gety();
					cout << x << " "<<y;
				}
				cout<<endl;
				cout<<v.size()<<endl;
			}
				
}
//Add coordinates to the correspondant leaf of the map
  bool bintree::add(float* coord, float* indexLeaf){
	  int dim = -1;
	  do{
		  dim = (dim + 1) % 2;
		  
		  if(coord[dim]<=(this->root)){
			  
			  if(left==right)	break;
			  (this->left)->add(coord, &(this->root));
			  
		  }else{
			  
			  if(left==right)	break;
			  (this->right)->add(coord, &(this->root));
		  }
	  }
	  while (true);
	  
	  map<float*, leaf*>::iterator found = mapLeaves.find(indexLeaf);
	  
	  if(found==mapLeaves.end()){
		  leaf* l = new leaf ;
		  l->addToLeaf(coord);
		   mapLeaves.insert ( pair<float*,leaf*>(indexLeaf, l ));
	  }else{
		  ((found->second))->addToLeaf(coord);
	  }
	  return true;
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
bintree::create2DBST(float **array, int depth, int nbreserves, int left, int right)
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
	  t = new bintree(array[left][n]);
	  t->add(array[left], &(t->root));
	  return t; 
  }
  
  int med = (left + right) / 2;
  t->root = array[med][n];
  cout<<"root : "<<t->root<<endl;
  tl = create2DBST(array, depth, med, left, med);
  tr = create2DBST(array, depth, right, med + 1, right);
  t->left = tl;
  t->right = tr;
  
return t;
}

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
	  t = new bintree(array[left][n]);
	  t->add(array[left], &(t->root));
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
	  t = new bintree(array[left][n]);
	  t->add(array[left], &(t->root));
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
	
    cout<< "adresse tree dans niveles : "<<&tree<<endl;
	readTreeLeaves();
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
float dist (vector<float>& d1, vector<float>& d2){
    unsigned int a;
    float d,sum=0;
    for(a=0; a < d1.size(); a++){
      d=d2[a] - d1[a];
      sum += d*d;
    }
    return sqrt(sum);
}

float* bintree::searchNode(vector<float> bcoordarray, bintree* tree, int depth){
	unsigned int n = 0;
    depth++;
    if(depth%2 != 0) {n = 1;}
    cout<<(tree->root)<<" "<<(tree->left)<<" "<<(tree->right)<<endl;
	if((bcoordarray[n]<=(tree->root)) && (tree->left)!=(tree->right)){
		
		return searchNode(bcoordarray, tree->left, depth);
		
	}else if((tree->left)!=(tree->right)){
		
		return searchNode(bcoordarray, tree->right, depth);
	}
		return &(tree->root);
}

vector<float> bintree::evaluateDist(vector<float> bcoordarray, bintree* tree, vector<reserve> reserves){

    vector<float> rcoordarray;
    vector<float> distances;
     for (vector<reserve>::iterator k = reserves.begin(); 
                             k != reserves.end(); 
                             ++k) 
                {
					rcoordarray[0]=k->getx();
					rcoordarray[1]=k->gety();
				}
	distances.push_back(dist(rcoordarray, bcoordarray));
    return distances;
}

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

vector<float> findNearestReserve(int index, vector<reserve>& reserves){
    
    vector<float> output;
    
    vector<reserve>::iterator k = reserves.begin();
    advance (k,index); 
    output[0] = k->getx();
    output[1] = k->gety();
    return output;
}
vector<reserve> searchVectorNode(float* indexMap){
	vector<reserve> rcoordarray;
	for(auto it = mapLeaves.cbegin(); it != mapLeaves.cend(); ++it)
			{
				if((it->first)==indexMap){
					rcoordarray =(it->second)->getCoordLeaf();
					break;
				}
			}
	return rcoordarray;
}

void writeOutputFile(vector<float> output, ofstream &file){
	unsigned int a;
    for(a=0; a < output.size(); a++){
      file << output[a] << " ";
    }
    file << endl;
}

void bintree::readBases(ofstream &ofile, bintree* tree, string baseFile){
	
    cout<< "adresse tree dans readBases : "<<&tree<<endl;
    string bline;
    vector<float> output, cd, c;
    ifstream bfile ;
    bfile.open(baseFile);
    bool seen =false;
    unsigned int dim =2;
    
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
      float* indexMap = bintree::searchNode(c, tree, 0);
      //recuperation du vecteur de coordonnées de réserves stockées sur la feuille
      vector<reserve> reserves = searchVectorNode(indexMap);
      cd = bintree::evaluateDist(c, tree, reserves);
      output = findNearestReserve(min(cd), reserves);
      writeOutputFile(output, ofile);
      }

      bfile.close();
    }
    else cout << "Error when opening bases files " << endl;; 
  }

void bintree::openOutputFile(bintree* tree, string outputFile, string baseFile){
    
    cout<< "adresse tree dans opent output : "<<&tree<<endl;
    ofstream ofile;
    ofile.open (outputFile.c_str());
    if (!ofile.is_open())
    {
  	  cout << "Error when opening output file " << endl;;
    }
    readBases(ofile, tree, baseFile);
    ofile.close();
  }

  int main (int argc, char* argv[]) {

  	bool points=false, input=false, output=false, heuristic=false;
  	vector<reserve> reserves;
    float **reservesArray=NULL;
  	string basef, outputf, splitting;
    vector <string> cinbases;
    vector <float> outputv;
    bintree *tree ;
    cout<< "adresse tree debut main : "<<&tree<<endl;
    int depth = -1;
  	int nbreserves;
  	for (int i = 1; i < argc; ++i) {
          string arg = argv[i];
  		if((arg == "-p")||(arg=="--points")) {
  			points=true;
  			if (i + 1 < argc) {
  				nbreserves = readReserves(reserves, argv[++i]);
				reservesArray = setReservesArray(reserves);
				read(reservesArray);

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
    if(!points) {cerr << "ERROR : -p or --points is a compulsory option." << endl; exit (EXIT_FAILURE);}

    // if((!input) && (output))
  	// {
   //    cinbases = readBasesFromInput();
   //    ofstream myfile;
   //    myfile.open (outputf);
   //    readBases(reserves, outputf, cinbases, myfile);
   //    myfile.close();
   //    cout<<"The answers are in the output file you provided."<<endl;
    
  	// }else 
  	// if((input) && (output))
  	// {
  	//  openOutputFile(reserves,tree,basef);
  	// }//else if((input) && (!output))
  	// { 
   //   readBasesNoOutput(reserves, basef);
  	// }else if((!input) && (!output))
  	// { 
   //    cinbases = readBasesFromInput();
   //    readBasesNoInputOutput(reserves, cinbases);
	  // }
	if(heuristic){
		if(splitting == "mediana"){
			
			tree = bintree::create2DBST(reservesArray, depth, nbreserves);
			
			bintree::niveles(tree); //Recorrido: por niveles
			
			if (tree != nullptr) {delete tree;} //Destrucción
		}else if(splitting == "mitad"){
			bintree *tree = bintree::createhalf2DBST(reservesArray, depth, nbreserves);
			bintree::niveles(tree); //Recorrido: por niveles
			if (tree != nullptr) {delete tree;} //Destrucción
		}else if(splitting == "promedio"){
			bintree *tree = bintree::createfraction2DBST(reservesArray, depth, nbreserves);
			bintree::niveles(tree); //Recorrido: por niveles
			if (tree != nullptr) {delete tree;} //Destrucción
		}else{
			cerr << "-s option accepts only mediana, mitad or promedio arguments." << endl;
            exit (EXIT_FAILURE);
		}
	}else{
		
    cout<< "adresse tree avant create : "<<&tree<<endl;
		tree = bintree::create2DBST(reservesArray, depth, nbreserves);
    cout<< "adresse tree avant niveles : "<<&tree<<endl;
		bintree::niveles(tree); //Recorrido: por niveles
		if (tree != nullptr) {delete tree;} //Destrucción
	}
	if((input) && (output))
  	 {
		 
    cout<< "adresse tree avant outputfile : "<<&tree<<endl;
  	  bintree::openOutputFile(tree, outputf, basef);
  	 }
    for (unsigned int i = 0;i<reserves.size();i++) { delete[] reservesArray[i]; }
    delete[] reservesArray;
    
    return 0;
  }
