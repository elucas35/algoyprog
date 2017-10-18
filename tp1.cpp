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
  static bintree *create2DBST(float** array, int depth, int left = 0, int right = -1);
  
};


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
    float ** reservesArray;
    reservesArray = new float *[reserves.size()];
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
    for(i=0; i<10; i++){
      for(j=0;j<(sizeof(*array)/sizeof(**array)) ;j++){
        cout << array[i][j] << " ";
      }
      cout << endl;
    }
    int r = sizeof array / sizeof *array;
    cout << "row : " << r << endl;
    cout << "col : " << sizeof *array / sizeof **array  << endl;

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

bintree *
bintree::create2DBST(float **array, int depth, int left, int right)
{
  int n=0;
  depth++;
  bintree *t = new bintree;
  bintree *tl;
  bintree *tr;

  if(depth%2 != 0) { n = 1; }
  
  //quickSort2D(array, left, right -1, n);
  if (right == -1) {right = 10;}
  if (left == right) { return nullptr; }
  if( left == right - 1 ) { return new bintree(array[left][n]); }
  

  int med = (left + right) / 2;
  t->root = array[med][n];
  tl = create2DBST(array, depth, left, med);
  tr = create2DBST(array, depth, med + 1, right);
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
  //quickSort2D(array, left, right -1, n);

  //if (left == right - 1 ) {cout << "med to one element is element " << v[left] << endl;}
      
      if(left != right - 1 && right != left){ // && (left!= right -2))
          cout << "sorted coordinates in "<< n << endl;
          for(i=left; i<right; i++){
            for(j=0;j<2;j++){
              cout << array[i][j] << " ";
            }
            cout << endl;
          }
      int med = (left + right) / 2;
      cout << "median is" << array[med][n]  << " and depth " << d <<endl;
      cout << "left is : " << left << "and ";
      cout << "right is : " << right << endl << endl;
      findmedian2D(left, med, d, array);
      findmedian2D(med +1, right, d, array);
    }
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

  	bool points=false, input=false, output=false;
  	vector<reserve> reserves;
    float **reservesArray;
  	string basef, outputf, line;
    vector <string> cinbases;
    vector <float> outputv;
    float ints[] = {16,2,77,29};
    int depth = -1;
  	
  	for (int i = 1; i < argc; ++i) {
          string arg = argv[i];
  		if((arg == "-p")||(arg=="--points")) {
  			points=true;
  			if (i + 1 < argc) {
  				readReserves(reserves, argv[++i]);
          reservesArray = setReservesArray(reserves);
          //read(reservesArray);

          bintree *tree = bintree::create2DBST(reservesArray, depth);
          bintree::niveles(tree); //Recorrido: por niveles
          if (tree != nullptr) {delete tree;} //Destrucción

          //findmedian2D(0,9, depth, reservesArray);

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

    for (int i = 0;i<reserves.size();i++) { delete[] reservesArray[i]; }
    delete[] reservesArray;
    return 0;
  }
