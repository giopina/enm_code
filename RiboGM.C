/* Program RiboGM:                  */
/* compute a gaussian network model */
/* for a given RNA molecule.        */
/* Created by Giovanni Pinamonti    */
/* PhD student at SISSA, Trieste    */
/* November 20th, 2013              */


#include <iostream>

#include "Structure3d.h"
#include "ElasticNet.h"

using namespace std;


int main(int argc, char*argv[]){

  if (argc<1) {cout<<"WARNING: not enough input parameters"<<endl; return 0;}
  
  string name(argv[1]);

  ElasticNet ENM;
  cout<<"Reading the parameters..."<<endl;
  ENM.readParameters("PARAMS_RIBOGM.DAT");


  cout<<"Reading the pdb file..."<<endl;
  ENM.readPDBFile(name);
  //  Structure3d structure;
  //  vector<string> list;
  //  list.push_back("P");
  //  list.push_back("C2");
  //  list.push_back("C4'");

  //  structure.setBeadList(list);
  //  structure.readFromPDBFile(argv[1]);
  cout<<"Constructing the gaussian model..."<<endl;
  //  ENM.setStructure(structure);
  
  cout<<"...the contact map..."<<endl;
  ENM.constructContactMap();
  cout<<"...and the interaction matrix!"<<endl<<endl;
  ENM.constructIntMat();
  
  cout<<"Let's go!"<<endl;
  ENM.Solve();
  cout<<"...done!"<<endl<<endl;
  
  cout<<"OUTPUT:"<<endl;
  ENM.Dump(name);
  //ENM.dumpMSD(name);
  
  cout<<"Done!"<<endl;
  
  ENM.computeForCovMatrix();
  ENM.dumpForCovMatrix("prova_for");
  
  return 0;
}
