#ifndef INPUTOUTPUT_H_INCLUDED
#define INPUTOUTPUT_H_INCLUDED

#include "bintree.h"
using namespace std;

void writeOutputFile(vector<float> output, ofstream &file);
void writeOutput(vector<float> output);
vector<float> readCoord(string line, bool &seen);
vector<float> readCoord(string line, bool &seen, string &filename);
void readBasesNoOutput(bintree* tree, string baseFile, int nbreserves);
vector<string> readBasesFromInput();
void readBases(bintree* tree, string outputFile, vector<string>& cinbases, ofstream& myfile, int nbreserves);
void readBasesNoInputOutput(bintree* tree, vector<string>& cinbases, int nbreserves);
void readBases(ofstream &ofile, bintree* tree, string baseFile, int nbreserves);
void openOutputFile(bintree* tree, string outputFile, string baseFile, int nbreserves);
#endif
