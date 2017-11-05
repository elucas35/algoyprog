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


//map (leaf's id, leaf of the tree)
map<bintree*,leaf*> mapLeaves ;


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
				read(reservesArray, nbreserves);
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
			
		}else if(splitting == "mitad"){
			bintree *tree = bintree::createhalf2DBST(reservesArray, depth, nbreserves);
			bintree::niveles(tree); //Recorrido: por niveles
			
		}else if(splitting == "promedio"){
			bintree *tree = bintree::createfraction2DBST(reservesArray, depth, nbreserves);
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
  	 if (tree != nullptr) {delete tree;} //Destrucción
    for (unsigned int i = 0;i<reserves.size();i++) { delete[] reservesArray[i]; }
    delete[] reservesArray;
    
    return 0;
}
