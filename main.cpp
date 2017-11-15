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
#include "reserve.h"
#include "bintree.h"
#include "inputoutput.h"
using namespace std;


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





//adds a leaves
void addLeavesToTree(float** reservesArray, bintree* t, int nbreserves){
  int i=0;
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
  int i=0;
  unsigned j=0;
  for(i=0; i<nbreserves; i++){
    for(j=0;j<(sizeof(*array)/sizeof(**array)) ;j++){
      cout << array[i][j] << " ";
    }
    cout << endl;
  }
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
      cout<<"dimension "<<dimension<<endl;
      while(getline(file, line)){

        nbreserves++;
        c = readCoord(line, seen, reserveFile);
        if(c.size() != dimension) {
			cout<<"c "<<c.size()<<endl;
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

   for (int i = 1; i < argc; ++i) {
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

  if(!points) {cerr << "ERROR : -p or --points is a compulsory option." << endl; exit (EXIT_FAILURE);}
  if(heuristic){
    if(splitting == "mediana"){
     tree = bintree::create2DBST(reservesArray, depth, nbreserves);
   //  bintree::checkTree(tree);
     addLeavesToTree(reservesArray, tree, nbreserves);
     attributeLimitLeaves();

   }else if(splitting == "mitad"){
     tree = bintree::createhalf2DBST(reservesArray, depth, nbreserves);
    // bintree::checkTree(tree);
     addLeavesToTree(reservesArray, tree, nbreserves);
     attributeLimitLeaves();
	//bintree::readLeaves();
   }else if(splitting == "promedio"){
     tree = bintree::createfraction2DBST(reservesArray, depth, nbreserves);
    // bintree::checkTree(tree);
     addLeavesToTree(reservesArray, tree, nbreserves);
     attributeLimitLeaves();
    // bintree::readLeaves();
   }else{
     cerr << "-s option accepts only mediana, mitad or promedio arguments." << endl;
     exit (EXIT_FAILURE);
   }
  }else{

    tree = bintree::create2DBST(reservesArray, depth, nbreserves);
    addLeavesToTree(reservesArray, tree, nbreserves);
    attributeLimitLeaves();
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
  if (tree != nullptr) {delete tree;} 
  for (unsigned int i = 0;i<reserves.size();i++) { delete[] reservesArray[i]; }
    delete[] reservesArray;

return 0;
}
