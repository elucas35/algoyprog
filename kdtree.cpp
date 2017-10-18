#include <vector>
#include <cstdlib>
#include <iostream>
#include "queue.h" 
#include <fstream>
#include <string>
#include <iterator>
#include <math.h>
#include <cstring>
using namespace std;

  class reserve
  {
    private:
    vector<float> coordinates;
    
    public:
    void addACoordinate(float);
    vector<float> getCoordinates();
    int getDim();
  };

  void reserve::addACoordinate(float c){
    coordinates.push_back(c);
  }
  vector<float> reserve::getCoordinates(){
    return coordinates;
  }
  int reserve::getDim(){
    return coordinates.size();
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
          if(isdigit(cstr[i]) || cstr[i] == '.') s += cstr[i];
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
    unsigned int j, m;
      for (m=0; m<coordarray.size(); m++){
      r.addACoordinate(coordarray[m]);

      }
    reserves.push_back(r);
  }
  void readReserves(vector<reserve>& reserves, string reserveFile){
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
  }

class bintree {
private:
  float root;
  bintree *left;
  bintree *right;

public:
  bintree() : left(nullptr), right(nullptr) {
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
  static bintree *dameBST( vector <float> &vx,  vector <float> &vy, int depth, int i = 0, int j = -1);
};

void passFunc(float **a)
{
    int i, j;
    for(i=0; i<10; i++){
      for(j=0;j<2;j++){
        a[i][j] = rand() % 100;
      }
    }
}
// cout << array[6][0] << endl; // output ligne 6 col 0 
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
void quickSort(vector <float>& arr, float left, float right) {
      
      float i = left, j = right;
      float tmp;
      float pivot = arr[(left + right) / 2];
      //cout << "median is" << pivot  << endl;

 
      /* partition */
      while (i <= j) {
            while (arr[i] < pivot)
                  i++;
            while (arr[j] > pivot)
                  j--;
            if (i <= j) {
                  tmp = arr[i];
                  arr[i] = arr[j];
                  arr[j] = tmp;
                  i++;
                  j--;
            }
      };
 
      /* recursion */
      if (left < j)
            quickSort(arr, left, j);
      if (i < right)
            quickSort(arr, i, right);
}

bintree *
bintree::dameBST( vector <float> &vx,  vector <float> &vy, int depth, int i, int j)
{
  //quickSort(vx, 0, vx.size() -1);
  //quickSort(vy, 0, vy.size() -1);
  depth++;
  //if (depth == 1 ) {vx.swap(vy);}
  if (j == -1)
    j = vx.size();
  if (i == j)
    return nullptr;
  if( i == j - 1 )
    return new bintree(vx[i]);

  bintree *t = new bintree;
  bintree *tl;
  bintree *tr;
  int medio = (i + j) / 2;
  
  t->root = vx[medio];
  tl = dameBST(vx, vy, depth, i, medio);
  tr = dameBST(vx, vy, depth, medio + 1, j);
  t->left = tl;
  t->right = tr;
  
  return t;
  
}

void findmedian2D(int left, int right, int d, float **array)
{
  
  int n = 0, i, j;
  d++;
 if(d%2 != 0) {n = 1;}
   
  //if (right == -1) { right = v.size();}
  quickSort2D(array, left, right -1, n);

  //if (left == right - 1 ) {cout << "med to one element is element " << v[left] << endl;}
      
      if(left != right - 1 && right != left){ // && (left!= right -2))
          cout << "soorted coordinates in "<< n << endl;
          for(i=left; i<right; i++){
            for(j=0;j<2;j++){
              cout << array[i][j] << " ";
            }
            cout << endl;
          }
      // cout << array[6][0] << endl; // output ligne 6 col 0 
      int med = (left + right) / 2;
      cout << "median is" << array[med][n]  << " and depth " << d <<endl;
      cout << "left is : " << left << "and ";
      cout << "right is : " << right << endl << endl;
      findmedian2D(left, med, d, array);
      findmedian2D(med +1, right, d, array);
    }
}
void findmedian(int left, int right, int d, vector < vector<float> > &rcoord)
{
  
  int n = 0, a;
  vector <float> v;
  d++;
 if(d%2 != 0) {n = 1;}
    cout << "coordinates x are : ";
        for (vector< vector<float> >::iterator l = rcoord.begin(); 
                             l != rcoord.end(); 
                             ++l) 
                {
                  
                  cout << l->at(0) << " ";

                }
                cout << endl << "coordinates y are : ";
        for (vector< vector<float> >::iterator l = rcoord.begin(); 
                             l != rcoord.end(); 
                             ++l) 
                {
                  
                  cout << l->at(1) << " ";

                }
                cout << endl << "coordinates v are : ";
    for (vector< vector<float> >::iterator l = rcoord.begin(); 
                             l != rcoord.end(); 
                             ++l) 
                {
                  v.push_back(l->at(n));
                  cout << l->at(n) << " ";

                }
                cout << endl;
  if (right == -1) { right = v.size();}
  quickSort(v, left, right -1);

  //if (left == right - 1 ) {cout << "med to one element is element " << v[left] << endl;}
       cout << "coordinates sorted v are : ";
      if(left != right - 1 && right != left){ // && (left!= right -2))
          
          for(a=left; a<right; a++){
              cout << v[a] << " ";
          }
          cout << endl;
      
      int med = (left + right) / 2;
      cout << "median is" << v[med]  << " and depth " << d <<endl;
      cout << "left is : " << left << "and ";
      cout << "right is : " << right << endl;
      findmedian(left, med, d, rcoord);
      findmedian(med +1, right, d, rcoord);
    }
}


void
bintree::niveles(bintree *tree)
{
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
}



  int main (int argc, char* argv[]) {

    bool points=false, input=false, output=false;
    vector<reserve> reserves;
    string basef, outputf, line;
    vector <string> cinbases;
    vector <float> outputv;
    int a, n , m;
    
    for (int i = 1; i < argc; ++i) {
          string arg = argv[i];
      if((arg == "-p")||(arg=="--points")) {
        points=true;
        if (i + 1 < argc) {
          readReserves(reserves, argv[++i]);
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
      }
    }
    if(!points) {cerr << "ERROR : -p or --points is a compulsory option." << endl; exit (EXIT_FAILURE);}

    //vector<float> rcoord;
    vector< vector<float> > rcoord;
    for (vector<reserve>::iterator k = reserves.begin(); 
                             k != reserves.end(); 
                             ++k) 
                {
                  rcoord.push_back(k->getCoordinates()) ;
                  
                }


//cout << " x coord" ;
vector< vector<float> >::iterator l = rcoord.begin() + 2;
                //cout << l->at(0) << endl;
l = rcoord.begin() + ((3 + 5) /2) ;
                //cout << l->at(0) << endl;


  // Construcción del árbol
  vector<float> vx {1, 2, 3, 4, 5};
  vector<float> vy {6, 7, 8, 9, 10};
  n=0;
 l = rcoord.begin() + ((0 + (rcoord.size() -1)) /2) ;
  int pivot = l->at(n);
  //cout << "pivot" << pivot << endl;
  int depth = -1;
  //quickSort(vx, 0, vx.size() -1);
    for (vector<vector<float>>::iterator g = rcoord.begin(); 
                             g != rcoord.end(); 
                             ++g) 
                {
                  //cout << g->at(0) << " ";
                  
                }

float **array;
array = new float *[10];
for(int i = 0; i <10; i++)
    array[i] = new float[2];

passFunc(array);
//cout << "row : " << (sizeof(array)/sizeof(*array)) << endl;
//cout << "col : " << (sizeof(*array)/sizeof(**array)) << endl;
   cout << "original coord : " << endl;
    int i, j;
    for(i=0; i<10; i++){
      for(j=0;j<2;j++){
        cout << array[i][j] << " ";
      }
      cout << endl;
    }
    cout << endl << endl;
// cout << array[6][0] << endl; // output ligne 6 col 0 


  findmedian2D(0,9, depth, array);
delete array;



  //bintree *tree = bintree::dameBST(vx, vy, depth);


  // Recorrido: por niveles
  //std::cout << "niveles: ";
  //bintree::niveles(tree);
  //std::cout << "\n";

  //Destrucción
  //if (tree != nullptr) {delete tree;}

  return 0;
}
