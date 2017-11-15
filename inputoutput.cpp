#include "inputoutput.h"
//writes the nearest reserve found into an output file
void writeOutputFile(vector<float> output, ofstream &file){
 unsigned int a=0;
 for(a=0; a < output.size(); a++){
  file << output[a] << " ";
}
file << endl;
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
//case 1 : call readBases() 
void openOutputFile(bintree* tree, string outputFile, string baseFile, int nbreserves){

  ofstream ofile;
  ofile.open (outputFile.c_str());
  if (!ofile.is_open())
  {
   cout << "Error when opening output file " << endl;;
 }
 readBases(ofile, tree, baseFile, nbreserves);
 ofile.close();
}
//case 1 : input file + output file 
void readBases(ofstream &ofile, bintree* tree, string baseFile, int nbreserves){

  string bline="";
  vector<float> output, cd, c, distToAreas;
  vector<int> nearestNeighbors;
  ifstream bfile ;
  bfile.open(baseFile);
  bool seen =false;
  unsigned int dim =2;
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
      
      //cout << "feuille initiale " << indext-> root << endl;
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
