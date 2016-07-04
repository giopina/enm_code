/*
 * NormalModes.cpp
 *
 *  Created on: Jul 5, 2012
 *      Author: gpolles
 */

#define MAX_PATH_LENGTH 1024

#include <iostream>
#include <cassert>
#include <fstream>
#include <string>
#include <sstream>
#include <stdexcept>
#include <cstdlib>
#include <cstring>
#include "NormalModes.h"

using namespace std;



NormalModes::NormalModes() : _numModes(0){
}

NormalModes::~NormalModes() {
}


void NormalModes::readFromFile(const std::string& basename, size_t numModes, size_t dim, bool printProgress = false) {
  char fname[MAX_PATH_LENGTH] ;
  sprintf(fname,"%s%s", basename.c_str(), "_eigenvalues.dat");
  ifstream eigvalFile(fname);
  if (!eigvalFile.is_open()){
    cerr << CLRLINE << "Fatal error. Could not open " << fname << endl;
    exit(1);
  }

  try{
    _eigenvalues.reserve(numModes);
  }catch(std::bad_alloc const&){
    cerr << CLRLINE << "NormalMode eigenvalues memory allocation failed!" << endl << flush;
    exit(1);
  }

  while(eigvalFile.good()){
    if (_eigenvalues.size() == numModes) break;
    int index;
    double val;
    eigvalFile >> index >> val;
    _eigenvalues.push_back(val);
  }

  if(_eigenvalues.size() < numModes){
    cerr << CLRLINE << "Error. Could read only " <<  _eigenvalues.size() << " eigenvalues. I wanted "<<numModes<< endl;
    exit(1);
  }

  try{
    _eigenvectors.resize(numModes);
    for (size_t i = 0; i < numModes; ++i) {
      
      sprintf(fname, "%s%s%ld%s",basename.c_str(),"_eigenvector_",i,".dat");
      ifstream eigvecFile(fname);
      if (!eigvecFile.is_open()){
        cerr << CLRLINE << "Fatal error. Could not open " << fname << endl<< flush;
        exit(1);
      }
      
      _eigenvectors[i].reserve(dim);
      for (size_t j = 0; j < dim; ++j) {
        int index;
        Vector3d p3d;
        eigvecFile >> index >> p3d;
        _eigenvectors[i].push_back(p3d);
	
      }
      if(eigvecFile.fail()){
        cerr << CLRLINE << "Fatal error reading " << fname << endl << flush;
        exit(1);
      }
      if(printProgress) printProgressBar(i,numModes,std::cout);
      
    }
    
  }catch(std::bad_alloc const&){
    cerr << CLRLINE << "NormalMode eigenvector memory allocation failed!" << endl << flush;
    exit(1);
  }


  _numModes = numModes;

}

real_t NormalModes::getEigenvalue(size_t i) const {
assert(i<_eigenvalues.size());
  return _eigenvalues[i];
}

const std::vector<Vector3d>& NormalModes::getEigenvector(size_t i) const {
  assert(i<_eigenvectors.size());
  return _eigenvectors[i];
}

void NormalModes::clear() {
  _eigenvalues.resize(0);
  _eigenvectors.resize(0);
  _numModes = 0;
}
