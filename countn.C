#include <iostream>
#include <string>
#include <cmath>
#include <boost/random.hpp>
#include "omp.h"
#include <vector>



#include "Structure3d.h"
#include "PrincipalComp.h"
#include "ElasticNet.h"
//#include "Matrix.h"
using namespace std;


void help_display(){
  cout<<"Help:"<<endl<<
    "-f --> input file  (required) format .pdb"<<endl<<
    "-par --> parameters file (DEFAULT = PARAMS_RIBOGM.DAT)"<<endl<<
    "-oname --> output name"<<endl<<
    "-h help"<<endl;;
  return;
}


int main(int argc, char* argv[]) {
  char oname[200];
  char file_name[200];
  char par_name[200];

  int i=1;
  int nopar=0;
  while(i<argc){
    if(!strncmp(argv[i],"-f",2)){
      i=i+1;
      if(i>argc-1){cout<<"ERROR: empty parameter"<<endl; help_display(); return 0;}
      sprintf(file_name,argv[i]);
      nopar+=1;
      cout<<"will generate a trajectory for the file "<<file_name<<endl;
    }
    else if(!strncmp(argv[i],"-par",4)){
      i=i+1;
      if(i>argc-1){cout<<"ERROR: empty parameter"<<endl; return 0;}
      sprintf(par_name,argv[i]);
      nopar+=1;
      cout<<"will read parameters from the file "<<par_name<<endl;
    }

    else if(!strncmp(argv[i],"-oname",5)){
      i=i+1;
      if(i>argc-1){cout<<"ERROR: empty parameter"<<endl; return 0;}
      sprintf(oname,argv[i]);
      nopar+=1;
      cout<<"\"My name is "<<oname<<".\""<<endl;
    }

    else if(!strncmp(argv[i],"-h",2)){
      help_display();
      return 0;
    }
    else{
      cout<<"WARNING: input n. "<<i<<" "<<argv[i]<<" ---> invalid parameter"<<endl;
    }
    i=i+1;
  }
  if(nopar<3){cout<<"ERROR: insufficient number of parameters"<<endl;
    help_display;
    return 0;
  }



  ElasticNet ENM;

    ENM.readParameters(par_name);
    int  error=ENM.readPDBFile(file_name);
    if(error!=0) {cout<<"ERROR READING THE FILE"<<endl;}// return 0;}

    cout<<"     contact map..."<<endl;
    ENM.constructContactMap();  
    
    ENM.dumpNearestNeighbours(oname);

    return 0;

}
