/* Program written by Giovanni Pinamonti
  PhD student at 
  Scuola Internazionale Superiori di Studi Avanzati, Trieste, Italy
   begin june 11th 2013 */

#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <malloc.h>
#include <string.h>   
#include <iostream>

#include "my_malloc.h" //libraries to allocate vectors and matrices
#include "io.h" //definitions of input and output subroutines
#include "CMatrix.h" //definition of the class CMatrix

using namespace std;

int main(int argc, char **argv){//PROGRAM MAIN
  char file_name[200];

  if(argc>3) {cout<<"### ERROR: too many arguments ###"<<endl; return 0;}
  
  sprintf(file_name,argv[1]);
  
  CMatrix Cij;
  
  Cij.TrajCompute(file_name);

  cout<<" Done!"<<endl;

  double trace=  Cij.GetTrace();
  
  cout<<"TRACCIA = "<<trace<<endl;

}//END PROGRAM MAIN
