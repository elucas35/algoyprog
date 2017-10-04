  #include <iostream>
  #include <fstream>
  #include <string>
  #include <vector>
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


  float dist (vector<float>& d1, vector<float>& d2){
    int a;
    float d,sum=0;
    for(a=0; a < d1.size(); a++){
      d=d2[a] - d1[a];
      sum += d*d;
    }
    return sqrt(sum);
  }
  vector<float> readCoord(string line, bool &seen){
	char *cstr = new char[line.length() + 1];
    int i, l;
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
                  cout << "ERROR : Bases coordonates must contain only the dot character which is '.' and/or numbers." << endl;
                  cout << "Please replace the character " << cstr[i] << "."  << endl;
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
    int i, l;
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
                  cout << "Please replace the character " << cstr[i] << "."  << endl;
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
  void setReserveCoord(vector<float> coordarray, vector<reserve>& reserves){

    reserve r;
    int j, m;
      for (m=0; m<coordarray.size(); m++){
      r.addACoordinate(coordarray[m]);

      }
      reserves.push_back(r);

     for(j=0;j<reserves.size();j++){
      reserves[j].getCoordinates();
      }
  }
  vector<float> evaluateDist(vector<float> bcoordarray, vector<reserve>& reserves){

    vector<float> rcoordarray;
    vector<float> distances;
    int a;

    for (vector<reserve>::iterator k = reserves.begin(); 
                             k != reserves.end(); 
                             ++k) 
                {
                    rcoordarray = k->getCoordinates();
                    distances.push_back(dist(rcoordarray, bcoordarray));
                    
                }

    return distances;
  }

  int min(vector<float> &d){
    int a, indice = 0;
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
    output = k->getCoordinates();
    
    return output;
  }

  void writeOutputFile(vector<float> output, ofstream &file){
    int a;
    
    for(a=0; a < output.size(); a++){
      file << output[a] << " ";
    }
    file << endl;
  }
  void writeOutput(vector<float> output){
    int a;
    for(a=0; a < output.size(); a++){
      cout << output[a] << " ";
    }
    cout << endl;
  }
  void readBases(ofstream &ofile, vector<reserve>& reserves, string baseFile){
    string bline;
    vector<float> output, cd, c;
    ifstream bfile ;
    bfile.open(baseFile);
    bool seen =false;
    if (bfile.is_open())
    {
      while(getline(bfile, bline)){
  		c = readCoord(bline, seen, baseFile);
  		 if(!c.size()) {
  		 ofile.close();
  		 ofile.open("output.txt", std::ofstream::out | std::ofstream::trunc);
  		 ofile << "ERROR : One line is empty, please fill the empty line and re-run the program." << endl;
  		 cout << "ERROR : the file must not contain any empty line, please fill the empty line and re-run the program." << endl; exit (EXIT_FAILURE);
      }
      cd = evaluateDist(c, reserves);
      output = findNearestReserve(min(cd), reserves);
      writeOutputFile(output, ofile);
      }

      bfile.close();
    }
    else cout << "Error when opening bases files " << endl;; 
  }
  void readBasesNoOutput(vector<reserve>& reserves, string baseFile){
    cout << "The nearest reserves to the bases are : " << endl;
    string bline;
    vector<float> output, cd, c;
    ifstream bfile ;
    bfile.open(baseFile);
    bool seen =false;
    if (bfile.is_open())
    {
      while(getline(bfile, bline)){
      c = readCoord(bline, seen, baseFile);
      cd = evaluateDist(c, reserves);
      output = findNearestReserve(min(cd), reserves);
      writeOutput(output);
      }

      bfile.close();
    }
    else cout << "Error when opening bases files " << endl;; 
  }

  void readBasesNoInputOutput(vector<reserve>& reserves, vector<string>& cinbases){
	cout << "The nearest reserves to the bases are : " << endl;
    string bline;
    vector<float> output, cd, c;
    bool seen =false;
    int a;
    if(cinbases.size()){
		for(a=0;a<cinbases.size();a++){
		  c = readCoord(cinbases[a], seen);
		  cd = evaluateDist(c, reserves);
		  output = findNearestReserve(min(cd), reserves);
		  writeOutput(output);
		}
	}
  }
  void readBases(vector<reserve>& reserves, string outputFile, vector<string>& cinbases, ofstream& myfile){
    string bline;
    vector<float> output, cd, c;
    ifstream bfile ;
    bfile.open(outputFile);
    bool seen =false;
    int a, l;

    if(cinbases.size()){
    for(a=0;a<cinbases.size();a++){
      outputFile = "standard input";
      c = readCoord(cinbases[a], seen, outputFile);
       if(!c.size()) {
       cout << "ERROR : the file must not contain any empty line, please fill the empty line and re-run the program." << endl; exit (EXIT_FAILURE);
       }
      cd = evaluateDist(c, reserves);
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
      cd = evaluateDist(c, reserves);
      output = findNearestReserve(min(cd), reserves);
      }

      bfile.close();
    }
    else cout << "Error when opening bases files " << endl;;
  }


  }

  void openOutputFile(vector<reserve>& reserves, string outputFile, string baseFile){
    
    ofstream ofile;
    ofile.open (outputFile.c_str());
    if (!ofile.is_open())
    {
  	  cout << "Error when opening output file " << endl;;
    }
    readBases(ofile, reserves, baseFile);
    ofile.close();
  }
  void readReserves(vector<reserve>& reserves, string reserveFile){
    string line, d;
    int dimension, nbreserves = 0;
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
      if(!c.size()) {
            cout << "ERROR : the file "<< reserveFile << " must not contain any empty line, please fill the empty line and re-run the program." << endl; exit (EXIT_FAILURE);
          }
      setReserveCoord(c, reserves);
      }
      file.close();
    }

    else cout << "Error when opening reserves files " << endl;;
  }
  int main (int argc, char* argv[]) {
  	if (argc < 3) {
  		cerr << "Error arguments missing " << endl;
  	}
  	bool points=false;
  	bool input=false;
  	bool output=false;
  	vector<reserve> reserves;
  	string basef;
  	string outputf;
    int nbases, a=0;
    string line;
    vector <string> cinbases;
    vector <float> outputv;
  	
  	for (int i = 1; i < argc; ++i) {
          string arg = argv[i];
  		if((arg == "-p")||(arg=="--points")) {
  			points=true;
  			if (i + 1 < argc) {
  				readReserves(reserves, argv[++i]);
  			} else {
  			  cerr << "-p option requires one argument." << endl;
  			}
  		} else if((arg == "-i")||(arg=="--input")) {
  			input=true;
  			if (i + 1 < argc) { 
  				basef=argv[++i];
  			}else {
  			  cerr << "-i option requires one argument." << endl;
  			}
  		} else if((arg=="-o")||(arg=="--output")){ 
  			output=true;
  			if (i + 1 < argc) { 
  				outputf=argv[++i];
  			}else {
  			  cerr << "-o option requires one argument." << endl;
  			}
  		}
  	}
  	if(!input && output)
  	{
      cout<<"Please enter the number of bases you want to install: "<<endl;
      cin >> nbases;
      cout<<"You chose to install " << nbases << "." <<endl;
      cout<<"Please enter bases coordinates, one per line, using same dimensions as for reserves: "<<endl;
      while(a<=nbases){
		  getline(cin, line);
		  if (line.size()) cinbases.push_back(line);
		  a++;
      }
      ofstream myfile;
      myfile.open (outputf);
      readBases(reserves, outputf, cinbases, myfile);
      myfile.close();
      cout<<"The answers are in the output file you provided."<<endl;
    
  	}else if((input) && (output))
  	{
  	 openOutputFile(reserves,outputf,basef);
  	}else if((input) && (!output))
  	{ 
     readBasesNoOutput(reserves, basef);
  	}else if((!input) && (!output))
  	{
		cout<<"Please enter the number of bases you want to install: "<<endl;
      cin >> nbases;
        cout<<"You chose to install " << nbases << "." <<endl;
      cout<<"Please enter bases coordinates, one per line, using same dimensions as for reserves: "<<endl;
      while(a<=nbases){
		  getline(cin, line);
		  if (line.size()) cinbases.push_back(line);
		  a++;
      }
      readBasesNoInputOutput(reserves, cinbases);
	}
  	
  	if(!points) cerr << "-p or --points is a compulsory option." << endl; exit (EXIT_FAILURE);
    return 0;
  }
