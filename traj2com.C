/* Program written by Giovanni Pinamonti
  PhD student at 
  Scuola Internazionale Superiori di Studi Avanzati, Trieste, Italy
   begin decembetr 10th 2013 */

#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <malloc.h>
#include <string.h>   
#include <iostream>


#include "Structure3d.h"
//#include "Matrix.h" //definition of the class CMatrix
//#include "PrincipalComp.h"

using namespace std;


int main(int argc, char **argv){//PROGRAM MAIN


  Structure3d structure(0);

  char input_name[200], output_name[200];
  sprintf(input_name,argv[1]);  
  sprintf(output_name,argv[2]);

  FILE* input;
  if ((input = fopen (input_name, "r")) == NULL) {
    fprintf (stderr,"Could not open file %s.\n.Exiting.\n", input_name);
    exit (1);    }

  ofstream output(output_name);


  int Nsteps=0;
  int error=structure.centersFromPDBFile(input);

  while(error==0){
    ++Nsteps;
    output<<"MODEL        "<<Nsteps<<endl;
    if(Nsteps%100==0)
      cout<<"Step "<<Nsteps<<endl;//" dumping"<<endl;    
    structure.dumpPDBFile(output);
    output<<"ENDMDL"<<endl;
    error=structure.centersFromPDBFile(input);
  }//end while



}//END MAIN
