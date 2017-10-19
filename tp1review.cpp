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
	 // int id[];

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
	  //void readTree(float **array, int depth,int left, int right);
	  bool add(float* coord, float index);
};

	map<float,leaf*> mapLeaves ;
	
	
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
    for(int i = 0; i <reserves.size(); i++)
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

    int i, j;
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
	   for(auto it = mapLeaves.cbegin(); it != mapLeaves.cend(); ++it)
			{
				cout << "leaf : "<< it->first << " vector : " ;
				vector<reserve> v =((it->second)->getCoordLeaf());
				for(auto it2=v.cbegin();it2!=v.cend();++it2){
					reserve r = *it2;
					float x = r.getx();
					float y = r.gety();
					cout << x << " "<<y;
				}
				cout<<endl;
			}
}
  bool bintree::add(float* coord, float indexLeaf){
	  int dim = -1;
	  do{
		  dim = (dim + 1) % 2;
		  if(coord[dim]<=(this->root)){
			  if(left==right)	break;
			  (this->left)->add(coord, root);
		  }else{
			  if(left==right)	break;
			  indexLeaf++;
			  (this->right)->add(coord, root);
		  }
	  }
	  while (true);
	  
	  map<float, leaf*>::iterator found = mapLeaves.find(indexLeaf);
	  
	  if(found==mapLeaves.end()){
		  leaf* l =new leaf ;
		  l->addToLeaf(coord);
		   mapLeaves.insert ( pair<float,leaf*>(indexLeaf, l ));
	  }else{
		  ((found->second))->addToLeaf(coord);
	  }
	 readTreeLeaves();
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
  
  int med = (left + right) / 2;
  
  if (left == right) { 
	  return nullptr; 
  }
  if( left == right - 1 ) {
	  t->add(array[left], array[left][n]);
	  return new bintree(array[left][n]); 
  }
  
  t->root = array[med][n];
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
	  t->add(array[left], array[left][n]);
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
	  t->add(array[left], array[left][n]);
	  return new bintree(array[left][n]);
  }
	  
  float fractionpoint;
  if(left != right - 1 && right != left){ 
	  cout << "sorted coordinates in level "<< n << endl;
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

  int main (int argc, char* argv[]) {

  	bool points=false, input=false, output=false, heuristic=false;
  	vector<reserve> reserves;
    float **reservesArray=NULL;
  	string basef, outputf, splitting;
    vector <string> cinbases;
    vector <float> outputv;
    
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
    
  	// }else if((input) && (output))
  	// {
  	//  openOutputFile(reserves,outputf,basef);
  	// }else if((input) && (!output))
  	// { 
   //   readBasesNoOutput(reserves, basef);
  	// }else if((!input) && (!output))
  	// { 
   //    cinbases = readBasesFromInput();
   //    readBasesNoInputOutput(reserves, cinbases);
	  // }
	if(heuristic){
		if(splitting == "mediana"){
			bintree *tree = bintree::create2DBST(reservesArray, depth, nbreserves);
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
		bintree *tree = bintree::create2DBST(reservesArray, depth, nbreserves);
		bintree::niveles(tree); //Recorrido: por niveles
		if (tree != nullptr) {delete tree;} //Destrucción
	}
	
    for (int i = 0;i<reserves.size();i++) { delete[] reservesArray[i]; }
    delete[] reservesArray;
    
    return 0;
  }
