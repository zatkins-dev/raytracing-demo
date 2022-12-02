//  ================================================================================
//  || ||
//  ||              ovenWalls ||
//  ||              ------------------------------------------------------ ||
//  ||              T H E R M A L   R A D I A T I O N ||
//  || ||
//  ||              D E M O N S T R A T I O N   C O D E ||
//  ||              ------------------------------------------------------ ||
//  || ||
//  ||       Developed by: Scott R. Runnels, Ph.D. ||
//  ||                     University of Colorado Boulder ||
//  || ||
//  ||                For: CU Boulder CSCI 4576/5576 and associated labs ||
//  || ||
//  ||           Copyright 2020 Scott Runnels ||
//  || ||
//  ||                     Not for distribution or use outside of the ||
//  ||                     this course. ||
//  || ||
//  ================================================================================

#include <cstdlib>
#include <ctime>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include "Cylinder.hpp"
#include "geom.hpp"

using std ::cout;
using std ::endl;
using std ::string;
using std ::stringstream;
using std ::vector;

typedef vector<double> VD;
typedef vector<vector<double>> VDD;
typedef vector<int> VI;
typedef vector<vector<int>> VII;

void FatalError(string msg) {

  cout << " " << endl;
  cout << " " << endl;
  cout << "Fatal Error: " << msg << endl;
  cout << " " << endl;
  exit(0);
}
